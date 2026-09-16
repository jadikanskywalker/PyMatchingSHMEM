#!/bin/bash
#SBATCH --job-name=pycall_gdb
#SBATCH --output=out/s-%j.out
#SBATCH --partition=zen4
#SBATCH --time=04:00:00

# gdb-wrapped variant of benchmark_shmem_obs_call.sh, for catching a stack trace on the mid-decode
# segfaults seen in the 36obs/72obs sweeps -- ASan is incompatible with Sandia OpenSHMEM's PGAS
# model, and PRTE's own --mca prte_stacktrace_output was confirmed not to fire for SOS binaries
# (see debug_tmp/scripts/debug_36obs_crash_native_diag.sh's own comment on this). Whole-process
# gdb wrap (pe_gdb_wrapper.sh) is what actually caught a real backtrace before (job 144398, double-
# free inside Mwpm::shatter_blossom_and_extract_matches) -- known risk: PRTE's own liveness
# detection can kill a rank as "hung" under ptrace overhead before it ever reaches a crash (seen
# once, on a specific heavy config) -- a job dying that way is a false negative, not a real crash.

if [ $# -le 8 ]
  then
    echo "Args: [ntasks] [sockets] [ntasks_per_node] [nthreads] [M] [k] [dem] [det] [flips]"
    exit 1
else
    ntasks=$1
    sockets=$2
    ntasks_per_node=$3
    nthreads=$4
    M=$5
    k=$6
    dem=$7
    det=$8
    flips=$9
fi

suffix=M${M}_ntasks${ntasks}_sockets${sockets}_ntps${ntasks_per_node}_nthreads${nthreads}_k${k}_${SLURM_JOB_ID}
out=out_$suffix

mkdir $out
cd $out
log=log_shmem.out
preds=../preds/preds_$suffix.01

echo "Job ID: $SLURM_JOB_ID" >> $log

source ~/.bash_profile
conda activate pymatching

export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=16G

export OMP_NUM_THREADS=$nthreads
export OMP_PLACES="cores($nthreads)"
export OMP_PROC_BIND=true

threads_per_task=$nthreads

# See benchmark_shmem_obs_call.sh's identical comment for why PE= has to be dropped (not just
# re-bound) for tasks spanning more than one 128-core socket.
map_by=ppr:$ntasks_per_node:node:PE=$threads_per_task
bind_to=core
if [ $threads_per_task -gt 128 ]; then
    map_by=ppr:$ntasks_per_node:node
    bind_to=none
fi

start_parallel=$(date +%s)
oshrun  \
    -n $ntasks \
    --map-by $map_by \
    --bind-to $bind_to \
    --report-bindings \
    ~/PyMatchingSHMEM/debug_tmp/scripts/pe_gdb_wrapper.sh \
    ~/PyMatchingSHMEM/build_sos/pymatching predict \
        --dem $dem \
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
