#!/bin/bash
#SBATCH --job-name=scorep_variance
#SBATCH --partition=zen4
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=64G

# Whole-node reservation (--exclusive) rather than --sockets=1 -- these thread counts (max 64) fit
# comfortably in one zen4 socket (128 cores) regardless, --exclusive just also rules out another
# job's memory-bandwidth/LLC contention on the same node without having to reason per-job about it.

if [ $# -lt 4 ]
  then
    echo "Args: [mode: trace|profile] [nthreads] [detection_events.b8] [error_model.dem]"
    exit 1
else
    mode=$1
    nthreads=$2
    det=$3
    dem=$4
fi

source ~/.bash_profile
conda activate pymatching

scorep_build=~/PyMatchingSHMEM/build_threads_scorep/pymatching

export OMP_NUM_THREADS=$nthreads
export OMP_PLACES=cores
export OMP_PROC_BIND=close

export SCOREP_EXPERIMENT_DIRECTORY=scorep_results
export SCOREP_ENABLE_PROFILING=true
# Default (~16MB) was too small for this workload's call-tree volume -- confirmed via "Error: No
# free memory page available... Please increase SCOREP_TOTAL_MEMORY" aborting every p=0.01 run and
# the p=0.001 trace run outright. 1GB comfortably covers 100-1000 shots x up to 64 locations with
# the new fine-grained leaf_decode/fusion_decode regions.
export SCOREP_TOTAL_MEMORY=1G
if [ "$mode" == "trace" ]; then
    export SCOREP_ENABLE_TRACING=true
else
    export SCOREP_ENABLE_TRACING=false
fi

$scorep_build predict \
    --dem "$dem" \
    --in "$det" \
    --in_format b8 \
    --out predicted_obs_flips.01 \
    --out_format 01 \
    --rounds_per_partition 32 \
    --use_threads \
    > run.out 2> run.err

echo "exit: $?"
echo "Score-P result: $(pwd)/scorep_results"
