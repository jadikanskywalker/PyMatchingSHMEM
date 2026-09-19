#!/bin/bash
#SBATCH --job-name=pycall
#SBATCH --output=out/s-%j.out
#SBATCH --partition=zen4
#SBATCH --time=02:00:00
# Known-bad nodes (comma-separated): rpc-96-1's ~/sw/$PLATFORM resolved to a dir without
# oshrun at job start (2026-09-16, jobs 155973/155978) -- PATH silently missing it isn't
# fatal to the script (no set -e), so the job "succeeds" with 0/N predictions. Add more
# node names here, comma-separated, as they're confirmed bad.
#SBATCH --exclude=rpc-96-1

if [ $# -le 8 ]
  then
    echo "Args: [ntasks] [sockets] [ntasks_per_node] [nthreads] [M] [k] [dem] [det] [flips] [L (optional, default 16)]"
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
    L=${10:-16}   # extraction_unit_size
fi

# nodes=$((sockets / 1))
# if [ $nodes -lt 1 ]
#   then
#     nodes=0
# fi
suffix=M${M}_ntasks${ntasks}_sockets${sockets}_ntps${ntasks_per_node}_nthreads${nthreads}_k${k}_${SLURM_JOB_ID}
out=out_$suffix

mkdir $out
cd $out
log=log_shmem.out
preds=../preds/preds_$suffix.01

if [ ! -d preds ]
  then
    mkdir preds
fi

echo "Job ID: $SLURM_JOB_ID" >> $log

source ~/.bash_profile
conda activate pymatching

export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=16G

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
oshrun  \
    -n $ntasks \
    --map-by $map_by \
    --bind-to $bind_to \
    --report-bindings \
    ~/PyMatchingSHMEM/scripts/pe_output_wrapper.sh \
    ~/PyMatchingSHMEM/build_sos/pymatching predict \
        --dem $dem \
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

# echo
# echo Shots with differring predictions:
# awk 'NR==FNR{a[NR]=$0; n=NR; next} {
#   if (FNR>n || $0!=a[FNR]) { print FNR-1; out=1 }
# } END {
#   if (n>FNR) { for (i=FNR+1;i<=n;i++) { print i-1; out=1 } }
#   if (!out) print "no differences"
# }' predicted_obs_flips__threads.01 predicted_obs_flips__shmem.01
