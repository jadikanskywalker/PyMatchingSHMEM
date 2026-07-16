#!/bin/bash
#SBATCH --job-name=pymatching_call
#SBATCH --partition=zen4
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --sockets=1
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=8G

thisM=$1
thisThreads=$2

source ~/.bash_profile
conda activate pymatching

threads_build=~/PyMatchingSHMEM/build_threads/pymatching
# sos_build=~/PyMatchingSHMEM/build_sos/pymatching
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=$thisThreads
# export SHMEM_SYMMETRIC_SIZE=8G

start=$(date +%s)
$threads_build predict \
    --dem ../error_model.dem \
    --in ../detection_events.b8 \
    --in_format b8 \
    --out ../predicted_obs_flips.01 \
    --out_format 01 \
    --rounds_per_partition $thisM \
    --use_threads \
    --num_repeats 10 \
    --extraction_unit_size 8

# $SWHOME/sos_1.5_scalable/bin/oshrun  \
#     -n 1 \
#     --map-by ppr:1:package:PE=$thisThreads \
#     --bind-to core \
#     --report-bindings \
# $sos_build predict \
#     --dem error_model.dem \
#     --in detection_events.b8 \
#     --in_format b8 \
#     --out predicted_obs_flips.01 \
#     --out_format 01 \
#     --rounds_per_partition $thisM \
#     --cross_rank_fusion_window_size 1 \
#     --use_threads \
#     --num_repeats 10
end=$(date +%s)
echo "M${thisM} ${thisThreads}threads: $((end - start)) seconds" >> bench.out
