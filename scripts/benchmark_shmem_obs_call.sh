#!/bin/bash
#SBATCH --job-name=pycall
#SBATCH --output=out/s-%j.out
#SBATCH --partition=zen4
#SBATCH --time=02:00:00

if [ $# -le 8 ]
  then
    echo "Args: [ntasks] [sockets] [ntasks_per_socket] [nthreads] [M] [k] [dem] [det] [flips]"
    exit 1
else
    ntasks=$1
    sockets=$2
    ntasks_per_socket=$3
    nthreads=$4
    M=$5
    k=$6
    dem=$7
    det=$8
    flips=$9
fi

# nodes=$((sockets / 1))
# if [ $nodes -lt 1 ]
#   then
#     nodes=0
# fi
suffix=M${M}_ntasks${ntasks}_sockets${sockets}_ntps${ntasks_per_socket}_nthreads${nthreads}_k${k}_${SLURM_JOB_ID}
log=log_$suffix.out
preds=preds/preds_$suffix.01

echo "Job ID: $SLURM_JOB_ID" >> $log

source ~/.bashrc
conda activate pymatching

export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=16G

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
    ../build_sos/pymatching predict \
        --dem $dem \
        --in $det \
        --in_format b8 \
        --out $preds \
        --out_format 01 \
        --rounds_per_partition $M \
        --obs_coors_included \
        --cross_rank_fusion_window_size $k \
        --task_division_strategy observable \
        --use_threads \
        --num_repeats 10 \
        &>> $log
end_parallel=$(date +%s)
parallel_time=$((end_parallel - start_parallel))
echo "  NTasks=$ntasks Sockets=$sockets NTPSocket=$ntasks_per_socket Threads=$nthreads: $parallel_time seconds" >> bench.out

python3 ../scripts/combine_results.py $preds $ntasks

echo >> $log
echo SHMEM  >> $log
echo correct predictions: >> $log
paste -d " " $preds $flips | grep "1 1\|0 0" | wc -l >> $log
echo wrong predictions: >> $log
paste -d " " $preds $flips | grep "0 1\|1 0" | wc -l >> $log

# echo
# echo Shots with differring predictions:
# awk 'NR==FNR{a[NR]=$0; n=NR; next} {
#   if (FNR>n || $0!=a[FNR]) { print FNR-1; out=1 }
# } END {
#   if (n>FNR) { for (i=FNR+1;i<=n;i++) { print i-1; out=1 } }
#   if (!out) print "no differences"
# }' predicted_obs_flips__threads.01 predicted_obs_flips__shmem.01
