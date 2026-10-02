#!/bin/bash
#SBATCH --job-name=pycall
#SBATCH --output=out/s-%j.out
#SBATCH --partition=zen4
#SBATCH --time=08:00:00
# Known-bad nodes (comma-separated): rpc-96-1's ~/sw/$PLATFORM resolved to a dir without
# rpc-91-7 (2026-09-25): UCX context creation stalls for minutes there (ibv_reg_mr / uct_ud_iface init), on both baseline and
# no-trim runs, while other nodes finish it quickly.
# oshrun at job start (2026-09-16, jobs 155973/155978) -- PATH silently missing it isn't
# fatal to the script (no set -e), so the job "succeeds" with 0/N predictions. Add more
# node names here, comma-separated, as they're confirmed bad.
#SBATCH --exclude=rpc-96-1,rpc-91-[1-20],rpc-94-3

# Graph-cache variant of benchmark_shmem_obs_call.sh -- see benchmark_shmem_obs_cached.sh for
# why (presets too large to keep a raw .dem around) and the M-must-match-the-cache caveat.

if [ $# -le 8 ]
  then
    echo "Args: [ntasks] [sockets] [ntasks_per_node] [nthreads] [M] [k] [graph_cache] [det] [flips] [L (optional, default 16)]"
    exit 1
else
    ntasks=$1
    sockets=$2
    ntasks_per_node=$3
    nthreads=$4
    M=$5
    k=$6
    graph_cache=$7
    det=$8
    flips=$9
    L=${10:-16}   # extraction_unit_size
fi

suffix=M${M}_ntasks${ntasks}_sockets${sockets}_ntps${ntasks_per_node}_nthreads${nthreads}_k${k}_${SLURM_JOB_ID}
out=out_$suffix

mkdir $out
cd $out
log=log_shmem.out
preds=../preds/preds_$suffix.01

if [ ! -d ../preds ]
  then
    mkdir ../preds
fi

echo "Job ID: $SLURM_JOB_ID" >> $log

source ~/.bash_profile
conda activate pymatching

# SHMEM_IMPL (exported by benchmark_shmem_obs_cached.sh): sos (default) or osss. Picks that implementation's launcher
# ($SOS_DIR / $OSSS_DIR from ~/.bash_profile) and the binary built against it (build_sos / build_osss).
shmem_impl=${SHMEM_IMPL:-sos}
case $shmem_impl in
    # impl_env: per-implementation environment, set on the command line through env so it reaches every PE whatever the
    # launcher forwards.
    #  SOS (libfabric verbs): SHMEM_OFI_PROVIDER=ofi_rxm and FI_VERBS_DEVICE_NAME=mlx5_2 are required (repacss
    #    test_sos_1.5.sh uses them; without the provider named, multi-PE init fails -- 2026-09-25: 8 and 16 real PEs on one
    #    node start with them). FI_MR_CACHE_MAX_COUNT=0 is SOS's own workaround for a libfabric memory-registration-cache
    #    bug with multithreaded SHMEM.
    #  OSSS (UCX): UCX_MEM_EVENTS=no -- without it 4+ PE runs fail with UCX memory-key deregistration errors
    #    (debug_tmp/osss_crash, 2026-09-25).
    sos)  shmem_dir=$SOS_DIR;  shmem_build=~/PyMatchingSHMEM/build_sos;  impl_env="SHMEM_OFI_PROVIDER=verbs;ofi_rxm FI_VERBS_DEVICE_NAME=mlx5_2 FI_MR_CACHE_MAX_COUNT=0" ;;
    osss) shmem_dir=$OSSS_DIR; shmem_build=~/PyMatchingSHMEM/build_osss; impl_env="UCX_MEM_EVENTS=no" ;;
    *)    echo "Unknown SHMEM_IMPL '$shmem_impl' (expected sos or osss)" >> $log; exit 1 ;;
esac
# SHMEM_BUILD_OVERRIDE (optional): use an alternate build dir (e.g. a RelWithDebInfo debug copy)
# instead of the impl's default build_sos/build_osss, without touching the default for other runs.
shmem_build=${SHMEM_BUILD_OVERRIDE:-$shmem_build}
# OSSS's oshrun is a python3.12 script; conda activate above replaces miniforge3/bin (which has it) on PATH.
export PATH="$shmem_dir/bin:$PATH:$HOME/miniforge3/bin"
# SURGERY_PRESET (exported by the driver): the seam spec the graph cache was generated from; it lets each PE build only the
# graph nodes it needs (the default subset build). Unset, the run falls back to a full graph build.
echo "SHMEM implementation: $shmem_impl ($shmem_dir), binary: $shmem_build/pymatching" >> $log

# Bigger than benchmark_shmem_obs_call.sh's 16G: these graphs (72/128/256/144obs) are larger
# than the 36obs one that value was tuned against -- regions_ptr/child_edges_ptr scale with
# the GLOBAL partition count (they do NOT shrink with more PEs, see gen_36obsx4.sh's own
# sizing notes), so this needs headroom regardless of ntasks/nthreads for this run.
#
# Confirmed real "Out of symmetric memory" failure (clean error, not a crash) at 32G for
# 256obs (18176 partitions, needed heap size 32G + ~25GB overrun before failing) -- 144obs has
# even more partitions (24108). Matching the largest gen_*obs.sh's own SHMEM_SYMMETRIC_SIZE
# (gen_144obs.sh uses 256G) rather than hand-tuning per preset again.
export SHMEM_SYMMETRIC_SIZE=48G

export OMP_NUM_THREADS=$nthreads
export OMP_PLACES="cores($nthreads)"
export OMP_PROC_BIND=true

threads_per_task=$nthreads

# Nodes here are dual-socket, 128 cores/socket. A task with more than 128 threads necessarily
# spans both sockets. map-by's PE=<N> clause itself mandates per-core binding (rejects any
# --bind-to but "core"/"hwt" -- "PE=<list> mapping directive cannot be combined with a binding
# directive other than core or hwt"), and --bind-to core alone rejects binding a single rank
# across packages ("failed to map... CPUs in more than one package") -- so PE= has to be dropped
# entirely for these, not just re-bound, leaving MPI's own placement unconstrained within the
# node and letting OMP_PLACES/OMP_PROC_BIND above do the actual in-process thread pinning.
map_by=ppr:$ntasks_per_node:node:PE=$threads_per_task
bind_to=core
if [ $threads_per_task -gt 128 ]; then
    map_by=ppr:$ntasks_per_node:node
    bind_to=none
fi

start_parallel=$(date +%s)
$shmem_dir/bin/oshrun  \
    -n $ntasks \
    --map-by $map_by \
    --bind-to $bind_to \
    --report-bindings \
    ~/PyMatchingSHMEM/scripts/pe_output_wrapper.sh \
    env $impl_env $shmem_build/pymatching predict \
        --graph_cache_path $graph_cache \
        --in $det \
        --in_format b8 \
        --out $preds \
        --out_format 01 \
        --rounds_per_partition $M \
        --cross_rank_fusion_window_size $k \
        --task_division_strategy observable \
        --extraction_unit_size $L \
        --extract_preemptively \
        --use_threads \
        ${SURGERY_PRESET:+--gen_surgery_preset $SURGERY_PRESET} \
        --num_repeats 10 \
    &>> $log

end_parallel=$(date +%s)
parallel_time=$((end_parallel - start_parallel))

echo "  NTasks=$ntasks Sockets=$sockets NTPNode=$ntasks_per_node Threads=$nthreads: $parallel_time seconds" >> ../bench.out

python3 ~/PyMatchingSHMEM/scripts/combine_results.py $preds $ntasks

echo >> $log
echo SHMEM  >> $log
# Multi-observable .01 format: one line per shot, num_observables bits concatenated per line --
# compare bit-by-bit (paste -d " " + grep only ever compares the two characters straddling the
# pasted space), mirroring scripts/compare_preds.sh's approach.
total_bits=$(awk '{ n += length($0) } END { print n }' $flips)
paste -d "" $preds $flips | awk -v total="$total_bits" '
{
    n = length($0) / 2
    for (i = 1; i <= n; i++) {
        if (substr($0, i, 1) == substr($0, n + i, 1)) correct++
        else wrong++
    }
}
END {
    print "correct predictions:"
    print correct + 0 "/" total
    print "wrong predictions:"
    print wrong + 0 "/" total
}' >> $log
