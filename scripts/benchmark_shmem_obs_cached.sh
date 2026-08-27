#!/bin/bash
#SBATCH --job-name=pybatch
#SBATCH --output=bench_obs-%j.out
#SBATCH --error=bench_obs-%j.err
#SBATCH --partition=zen4
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=1450GB

# Graph-cache variant of benchmark_shmem_obs.sh, for presets too large to keep a raw .dem
# around (dem_gen/scripts/gen_{72,128,256,144}obs.sh generate these via --graph_cache_path,
# never --dem_cache_path -- see project_dem_presets memory). Loads
# testdems/graph_<preset>_..._<rounds>r.cache directly via --graph_cache_path instead of
# --dem, skipping DEM parsing entirely.
#
# IMPORTANT: the graph cache's validity is keyed on enable_correlations, rounds_per_partition
# (M), and task_division_strategy (graph_cache.h:11-14) -- each cache was built with a specific
# M (task_division_strategy=observable, enable_correlations=false always; see the gen_*obs.sh
# scripts). Passing a mismatched M here will make read_user_graph_cache throw
# GraphCacheMismatchError, and unlike benchmark_shmem_obs.sh's DEM-based runs, there is no raw
# DEM to fall back to rebuilding from -- so M must match whichever cache [cache_suffix] selects
# (M21 for the plain caches from gen_{72,128,256,144}obs.sh, M22 for the _M22 caches from
# gen_{64,128,256}obs_m22.sh -- pass "_M22" as the 8th arg to select those). extraction_unit_size,
# extract_preemptively, cross_rank_fusion_window_size, and thread/PE counts are NOT part of
# the cache key and remain free to vary.

if [ $# -le 6 ]
  then
    echo "Args: [surgery_preset] [d] [p_dec] [shots] [rounds] [M] [k] [cache_suffix (optional)] [parent_dir (optional)]"
    exit 1
else
    surgery_preset=$1
    d=$2
    p_dec=$3
    shots=$4
    rounds=$5
    M=($6)
    k=$7
    cache_suffix=$8   # e.g. "_M22" -- gen_128obs_m22.sh/gen_256obs_m22.sh name their graph
                       # cache graph_<dem_suffix>_M22.cache (det/flips filenames are unaffected,
                       # only the graph cache filename carries this suffix) to distinguish the
                       # M22-aligned cache from the M21 one at a different round count sharing
                       # the same preset name; empty (default) matches the M21 caches' plain name.
fi

source ~/.bash_profile
conda activate pymatching

dem_suffix=${surgery_preset}_d${d}_p${p_dec}_${rounds}r
det_suffix=${dem_suffix}_${shots}s
# parent_dir (9th arg, optional): nests output under a named subfolder, e.g.
# bench_shmem_obs_many_observables, instead of directly under PyMatchingSHMEM/. The
# benchmark_shmem_obs_call_cached.sh path above is absolute specifically so this can nest to any
# depth without breaking the old "../scripts/..." assumption that dirname sits one level under
# PyMatchingSHMEM/.
parent_dir=${9:-.}
dirname=$parent_dir/bench_$det_suffix

mkdir -p $dirname

cd $dirname

if [ ! -d preds ]
  then
    mkdir preds
fi

if [ ! -d out ]
  then
    mkdir out
fi

graph_cache=~/PyMatchingSHMEM/testdems/graph_${dem_suffix}${cache_suffix}.cache
det=~/PyMatchingSHMEM/testdems/detection_events_$det_suffix.b8
flips=~/PyMatchingSHMEM/testdems/actual_obs_flips_$det_suffix.01
echo $graph_cache
echo $det
echo $flips


shmem_threads=(8 16 32 64 128)

# Fill sockets
shmem_ntasks2_sockets=(1 1 1 1 2)
shmem_ntasks2_ntps=(   0 0 2 2 1) # ntasks per socket

shmem_ntasks4_nodes=(  1 1 1 1 2)
shmem_ntasks4_sockets=(1 1 1 2 4)
shmem_ntasks4_ntps=(   0 4 4 2 1)

shmem_ntasks8_nodes=(  1 1 1 2 4)
shmem_ntasks8_sockets=(1 1 2 4 8)
shmem_ntasks8_ntps=(   8 8 4 2 1)

shmem_ntasks16_nodes=(   1 1 2 4 8)
shmem_ntasks16_sockets=( 1 2 4 8 16)
shmem_ntasks16_ntps=(   16 8 4 2 1)


shmem_quarter_threads=(16 32 32)
shmem_quarter_ntasks=( 4  4  2)
shmem_quarter_ntps=(   2  1  1)
shmem_quarter_sockets=(2  4  2)
shmem_quarter_nodes=(  1  2  1)


shmem_half_threads=(32 64 64)
shmem_half_ntasks=(  4  4  2)
shmem_half_ntps=(    2  1  1)
shmem_half_sockets=( 2  4  2)
shmem_half_nodes=(   1  2  1)


# Confirmed real OOM (MaxRSS ~268GB, killed and still climbing) at ntasks=1 with the old
# 256GB cap, for even the smallest-thread-count 72obs config -- graph load / partition-metadata
# construction alone need ~270GB+ per rank for these presets, independent of nthreads (all
# nthreads tested for a given preset/ntasks OOM'd identically). Since --exclusive already
# reserves the WHOLE node (~1.5TB here) regardless of --mem (no other job can share it), --mem
# only sets this job's own OOM-kill threshold -- there is no cost to just requesting close to
# the full node on every tier instead of hand-tuning a per-ntasks formula that keeps turning out
# to be wrong. NODE_MEM is per node (SLURM's --mem is per-node, not total, when --nodes>1), so
# this also correctly gives each node in a multi-node run its own ~1.45TB budget.
NODE_MEM=1450GB

repeats=1

serial_build=~/PyMatchingSHMEM/build/pymatching
shmem_build=~/PyMatchingSHMEM/build_sos/pymatching
echo "serial_build:  $serial_build" >> bench.out
echo "shmem_build: $shmem_build" >> bench.out
echo "----------" >> bench.out
echo "shots: $shots    rounds: $rounds" >> bench.out

for ((m=0; m<${#M[@]}; m++ )); do
    thisM=${M[$m]}
    echo "----------" >> bench.out
    echo "M: $thisM" >> bench.out
    echo "  SHMEM:" >> bench.out
    for ((i=0; i<${#shmem_threads[@]}; i++ )); do
        thisThreads=${shmem_threads[$i]}

        # # single socket runs
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=1 \
            --exclusive \
            --mem=$NODE_MEM \
            ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call_cached.sh \
                1 1 1 $thisThreads $thisM $k $graph_cache $det $flips
        done

        thisSockets=${shmem_ntasks2_sockets[$i]}
        thisNTPS=${shmem_ntasks2_ntps[$i]}
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=1 \
            --exclusive \
            --mem=$NODE_MEM \
            ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call_cached.sh \
                2 $thisSockets $thisNTPS $thisThreads $thisM $k $graph_cache $det $flips
        done

        thisNodes=${shmem_ntasks4_nodes[$i]}
        thisSockets=${shmem_ntasks4_sockets[$i]}
        thisSPN=$((thisSockets / thisNodes))
        thisNTPS=${shmem_ntasks4_ntps[$i]}
        thisNTPN=$((thisSPN * thisNTPS)) # ntasks per node
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=$thisNodes \
            --exclusive \
            --mem=$NODE_MEM \
            ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call_cached.sh \
                4 $thisSockets $thisNTPS $thisThreads $thisM $k $graph_cache $det $flips
        done

        # thisNodes=${shmem_ntasks8_nodes[$i]}
        # thisSockets=${shmem_ntasks8_sockets[$i]}
        # thisSPN=$((thisSockets / thisNodes))
        # thisNTPS=${shmem_ntasks8_ntps[$i]}
        # thisNTPN=$((thisSPN * thisNTPS)) # ntasks per node
        # thisMEM=$((8 * 160 / $thisNodes))
        # for ((r=0; r<repeats; r++)); do
        # sbatch \
        #     --nodes=$thisNodes \
        #     --exclusive \
        #     --mem=${thisMEM}GB \
        #     ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call_cached.sh \
        #         8 $thisSockets $thisNTPS $thisThreads $thisM $k $graph_cache $det $flips
        # done

        # thisNodes=${shmem_ntasks16_nodes[$i]}
        # thisSockets=${shmem_ntasks16_sockets[$i]}
        # thisSPN=$((thisSockets / thisNodes))
        # thisNTPS=${shmem_ntasks16_ntps[$i]}
        # thisNTPN=$((thisSPN * thisNTPS)) # ntasks per node
        # thisMEM=$((1280))
        # for ((r=0; r<repeats; r++)); do
        # sbatch \
        #     --nodes=$thisNodes \
        #     --exclusive \
        #     --mem=${thisMEM}GB \
        #     ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call_cached.sh \
        #         16 $thisSockets $thisNTPS $thisThreads $thisM $k $graph_cache $det $flips
        # done
    done
    for ((i=0; i<${#shmem_quarter_threads[@]}; i++ )); do
        # Quarter subscribe sockets
        thisThreads=${shmem_quarter_threads[$i]}
        thisNTasks=${shmem_quarter_ntasks[$i]}
        thisNodes=${shmem_quarter_nodes[$i]}
        thisSockets=${shmem_quarter_sockets[$i]}
        thisNTPS=${shmem_quarter_ntps[$i]}
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=$thisNodes \
            --exclusive \
            --mem=$NODE_MEM \
            ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call_cached.sh \
                $thisNTasks $thisSockets $thisNTPS $thisThreads $thisM $k $graph_cache $det $flips
        done

        # Half subsribe sockets
        thisThreads=${shmem_half_threads[$i]}
        thisNTasks=${shmem_half_ntasks[$i]}
        thisNodes=${shmem_half_nodes[$i]}
        thisSockets=${shmem_half_sockets[$i]}
        thisNTPS=${shmem_half_ntps[$i]}
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=$thisNodes \
            --exclusive \
            --mem=$NODE_MEM \
            ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call_cached.sh \
                $thisNTasks $thisSockets $thisNTPS $thisThreads $thisM $k $graph_cache $det $flips
        done
    done
done

start_serial=$(date +%s)
$serial_build predict \
    --graph_cache_path $graph_cache \
    --in $det\
    --in_format b8 \
    --out preds/preds_0.01 \
    --out_format 01 \
    --num_repeats 10 \
    >> log_0.out
end_serial=$(date +%s)
serial_time=$((end_serial - start_serial))
echo "0: $serial_time seconds" >> bench.out
