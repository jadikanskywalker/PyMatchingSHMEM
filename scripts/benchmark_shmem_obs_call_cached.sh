#!/bin/bash
#SBATCH --job-name=pycall
#SBATCH --output=out/s-%j.out
#SBATCH --partition=zen4
#SBATCH --time=08:00:00

# Graph-cache variant of benchmark_shmem_obs_call.sh -- see benchmark_shmem_obs_cached.sh for
# why (presets too large to keep a raw .dem around) and the M-must-match-the-cache caveat.

if [ $# -le 8 ]
  then
    echo "Args: [ntasks] [sockets] [ntasks_per_socket] [nthreads] [M] [k] [graph_cache] [det] [flips]"
    exit 1
else
    ntasks=$1
    sockets=$2
    ntasks_per_socket=$3
    nthreads=$4
    M=$5
    k=$6
    graph_cache=$7
    det=$8
    flips=$9
fi

suffix=M${M}_ntasks${ntasks}_sockets${sockets}_ntps${ntasks_per_socket}_nthreads${nthreads}_k${k}_${SLURM_JOB_ID}
out=out_$suffix

mkdir $out
cd $out
log=log_shmem.out
preds=../preds/preds_$suffix.01

echo "Job ID: $SLURM_JOB_ID" >> $log

source ~/.bash_profile
conda activate pymatching

export FI_VERBS_DEVICE_NAME="mlx5_2"
# Bigger than benchmark_shmem_obs_call.sh's 16G: these graphs (72/128/256/144obs) are larger
# than the 36obs one that value was tuned against -- regions_ptr/child_edges_ptr scale with
# the GLOBAL partition count (they do NOT shrink with more PEs, see gen_36obsx4.sh's own
# sizing notes), so this needs headroom regardless of ntasks/nthreads for this run.
#
# Confirmed real "Out of symmetric memory" failure (clean error, not a crash) at 32G for
# 256obs (18176 partitions, needed heap size 32G + ~25GB overrun before failing) -- 144obs has
# even more partitions (24108). Matching the largest gen_*obs.sh's own SHMEM_SYMMETRIC_SIZE
# (gen_144obs.sh uses 256G) rather than hand-tuning per preset again.
export SHMEM_SYMMETRIC_SIZE=256G

export OMP_NUM_THREADS=$nthreads
export OMP_PLACES="cores($nthreads)"
export OMP_PROC_BIND=true

threads_per_task=$nthreads

start_parallel=$(date +%s)
oshrun  \
    -n $ntasks \
    --map-by ppr:$ntasks_per_socket:package:PE=$threads_per_task \
    --bind-to core \
    --report-bindings \
    ~/PyMatchingSHMEM/scripts/pe_output_wrapper.sh \
    ~/PyMatchingSHMEM/build_sos/pymatching predict \
        --graph_cache_path $graph_cache \
        --in $det \
        --in_format b8 \
        --out $preds \
        --out_format 01 \
        --rounds_per_partition $M \
        --cross_rank_fusion_window_size $k \
        --task_division_strategy observable \
        --extraction_unit_size 4 \
        --extract_preemptively \
        --use_threads \
        --num_repeats 10 \
    &>> $log

end_parallel=$(date +%s)
parallel_time=$((end_parallel - start_parallel))

echo "  NTasks=$ntasks Sockets=$sockets NTPSocket=$ntasks_per_socket Threads=$nthreads: $parallel_time seconds" >> ../bench.out

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
