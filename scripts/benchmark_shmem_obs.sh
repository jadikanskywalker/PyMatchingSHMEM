#!/bin/bash
#SBATCH --job-name=pybatch
#SBATCH --output=bench_obs-%j.out
#SBATCH --error=bench_obs-%j.err
#SBATCH --partition=zen4
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=256GB

if [ $# -le 7 ]
  then
    echo "Args: [surgery_preset] [d] [p_dec] [shots] [rounds] [M] [k] [L] [parent_dir (optional)]"
    exit 1
else
    surgery_preset=$1
    d=$2
    p_dec=$3
    shots=$4
    rounds=$5
    M=($6)
    k=$7
    L=$8                  # extraction_unit_size, forwarded to benchmark_shmem_obs_call.sh
    parent_dir=${9:-.}   # e.g. bench_shmem_obs_many_observables, to keep output alongside the
                          # cached-graph presets' runs -- see benchmark_shmem_obs_cached.sh
fi

source ~/.bash_profile
conda activate pymatching

dem_suffix=${surgery_preset}_d${d}_p${p_dec}_${rounds}r
det_suffix=${dem_suffix}_${shots}s
dirname=$parent_dir/bench_${det_suffix}_L${L}

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

dem=~/PyMatchingSHMEM/testdems/error_model_$dem_suffix.dem
det=~/PyMatchingSHMEM/testdems/detection_events_$det_suffix.b8
flips=~/PyMatchingSHMEM/testdems/actual_obs_flips_$det_suffix.01
echo $dem
echo $det
echo $flips


shmem_threads=(16 32 64 128 256)

# Fill sockets. threads: 16 32 64 128 256 (index-aligned below). sockets is the TOTAL count of
# distinct sockets touched across the whole job (nodes * 2, since every node here is either fully
# packed onto one socket-worth of cores or spans both of its sockets -- see the 256-thread column,
# where a single ntasks_per_node=1 task at 256 threads necessarily spans both sockets of its node).
# Kept identical to benchmark_shmem_obs_cached.sh's own arrays -- see that script for the
# topology math these were validated against.
shmem_ntasks2_nodes=(  1 1 1 1 2)
shmem_ntasks2_sockets=(1 1 1 2 4)
shmem_ntasks2_ntpn=(   2 2 2 2 1) # ntasks per node

shmem_ntasks4_nodes=(  1 1 1 2 4)
shmem_ntasks4_sockets=(1 1 2 4 8)
shmem_ntasks4_ntpn=(   4 4 4 2 1)

shmem_ntasks8_nodes=(  1 1 2 4 8)
shmem_ntasks8_sockets=(1 2 4 8 16)
shmem_ntasks8_ntpn=(   8 8 4 2 1)


# NODE_MEM matches benchmark_shmem_obs_cached.sh's own flat-per-tier choice: --exclusive already
# reserves the whole node regardless of --mem, so there's no cost to requesting close to the full
# node on every tier instead of a fragile per-tier memory formula (even though these DEM-based
# presets need far less than the graph-cache ones that value was originally sized for).
NODE_MEM=1490GB

repeats=1

serial_build=~/PyMatchingSHMEM/build/pymatching
# threads_build=~/PyMatchingSHMEM/build_threads_release/pymatching
shmem_build=~/PyMatchingSHMEM/build_sos/pymatching
echo "serial_build:  $serial_build" >> bench.out
# echo "threads_build: $threads_build" >> bench.out
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
            ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call.sh \
                1 1 1 $thisThreads $thisM $k $dem $det $flips $L
        done

        thisNodes=${shmem_ntasks2_nodes[$i]}
        thisSockets=${shmem_ntasks2_sockets[$i]}
        thisNTPN=${shmem_ntasks2_ntpn[$i]}
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=$thisNodes \
            --exclusive \
            --mem=$NODE_MEM \
            ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call.sh \
                2 $thisSockets $thisNTPN $thisThreads $thisM $k $dem $det $flips $L
        done

        thisNodes=${shmem_ntasks4_nodes[$i]}
        thisSockets=${shmem_ntasks4_sockets[$i]}
        thisNTPN=${shmem_ntasks4_ntpn[$i]}
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=$thisNodes \
            --exclusive \
            --mem=$NODE_MEM \
            ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call.sh \
                4 $thisSockets $thisNTPN $thisThreads $thisM $k $dem $det $flips $L
        done

        thisNodes=${shmem_ntasks8_nodes[$i]}
        thisSockets=${shmem_ntasks8_sockets[$i]}
        thisNTPN=${shmem_ntasks8_ntpn[$i]}
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=$thisNodes \
            --exclusive \
            --mem=$NODE_MEM \
            ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call.sh \
                8 $thisSockets $thisNTPN $thisThreads $thisM $k $dem $det $flips $L
        done
    done
done

start_serial=$(date +%s)
$serial_build predict \
    --dem $dem \
    --in $det\
    --in_format b8 \
    --out preds/preds_0.01 \
    --out_format 01 \
    --num_repeats 10 \
    >> log_0.out
end_serial=$(date +%s)
serial_time=$((end_serial - start_serial))
echo "0: $serial_time seconds" >> bench.out