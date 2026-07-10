#!/bin/bash
#SBATCH --job-name=trace_obs
#SBATCH --output=trace-%j.out
#SBATCH --partition=zen4
#SBATCH --time=05:00:00
#SBATCH --mem=1000GB

echo $SLURM_JOB_ID

cd ~/PyMatchingSHMEM
source ~/.bash_profile
conda activate pymatching

export FI_VERBS_DEVICE_NAME="mlx5_2"

if [ $# -le 9 ]
  then
    echo "Args: [n] [pps] [package/node] [nthreads_shmem] [M] [k] [dir] [dem] [det] [flips]"
    exit 1
else
    n=$1
    pps=$2
    perwhat=$3
    nthreads_shmem=$4
    M=$5
    k=$6
    dir=$7
    dem=$8
    det=$9
    flips=${10}
fi

if [ ! -d $dir ]
  then
    mkdir $dir
fi

cd $dir

rm -r scorep_results* out_parallel 2>/dev/null
rm *.out *.err *.01 *.txt 2>/dev/null

export OMP_NUM_THREADS=$nthreads_shmem
export OMP_PLACES=cores
export OMP_PROC_BIND=true

# Profiling stays on for the cube_stat cross-check; tracing is the actual answer.
export SCOREP_ENABLE_PROFILING=true
export SCOREP_ENABLE_TRACING=true
export SCOREP_TOTAL_MEMORY=4G

echo "Starting SHMEM run..."
start_parallel=$(date +%s)
export SHMEM_SYMMETRIC_SIZE=16G
$SWHOME/sos_1.5_scalable/bin/oshrun  \
    -n $n \
    --map-by ppr:$pps:$perwhat:PE=$nthreads_shmem \
    --bind-to core \
    --report-bindings \
    ~/PyMatchingSHMEM/scripts/scorep_wrapper.sh \
    ~/PyMatchingSHMEM/build_sos_profile/pymatching predict \
    --dem ../testdems/$dem \
    --in ../testdems/$det \
    --in_format b8 \
    --out predicted_obs_flips__shmem.01 \
    --out_format 01 \
    --rounds_per_partition $M \
    --obs_coors_included \
    --cross_rank_fusion_window_size $k \
    --task_division_strategy observable \
    --use_threads \
    --num_repeats 10 \
    &> log_shmem.out

end_parallel=$(date +%s)
parallel_time=$((end_parallel - start_parallel))
echo "SHMEM run completed in $parallel_time seconds."

python3 ../scripts/combine_results.py predicted_obs_flips__shmem.01 $n

# Check work
echo SHMEM
echo correct predictions:
paste -d " " predicted_obs_flips__shmem.01 ../testdems/$flips | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips__shmem.01 ../testdems/$flips | grep "0 1\|1 0" | wc -l

echo "--- Wall-clock phase analysis ---"
conda run -n pymatching python3 ../scripts/analyze_trace_wall_clock.py scorep_results_pe*/traces.otf2 > analyze_trace_wall_clock.out 2>&1
tail -80 analyze_trace_wall_clock.out

cd ..
