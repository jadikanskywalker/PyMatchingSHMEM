#!/bin/bash
#SBATCH --job-name=PyMatchingSHMEM
#SBATCH --partition=zen4
#SBATCH --time=02:00:00

if [ $# -le 8 ]
  then
    echo "Args: [nodes] [processes_per_socket] [num_sockets_per_node] [nthreads_shmem] [M] [k] [dem] [det] [flips]"
    exit 1
else
    nodes=$1
    pps=$2
    num_sockets_per_node=$3
    nthreads=$4
    M=$5
    k=$6
    dem=$7
    det=$8
    flips=$9
fi

n=$((nodes * num_sockets_per_node * pps))
suffix=M${M}_n${n}_pps${pps}_${nthreads}_k${k}
log=log_$suffix.out
preds=preds_$suffix.01

source ~/.bashrc
conda activate pymatching

export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=12G

export OMP_NUM_THREADS=$nthreads

start_parallel=$(date +%s)
oshrun  \
    -n $n \
    --map-by ppr:$pps:package:PE=$nthreads \
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
        &> $log
end_parallel=$(date +%s)
parallel_time=$((end_parallel - start_parallel))
echo "  N=$n PPN=$pps (Nodes=$nodes) Threads=$nthreads: $parallel_time seconds" >> bench.out

python3 ../scripts/combine_results.py $preds $n

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
