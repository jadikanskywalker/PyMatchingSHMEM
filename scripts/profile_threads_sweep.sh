#!/bin/bash
#SBATCH --job-name=profile_threads_sweep
#SBATCH --output=profile_threads_sweep-%j.out
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G

if [ $# -lt 1 ]
  then
    echo "Args: [dir]"
    exit 1
else
    dir=$1
fi

source ~/.bash_profile
conda activate pymatching

mkdir -p $dir
cd $dir

d=21
rounds=10000
shots=1000
M=32
threads=(16 32 64 128 256)

stim gen \
    --rounds=$((rounds-1)) \
    --distance=$d \
    --after_clifford_depolarization=0.001 \
    --code surface_code \
    --task rotated_memory_x \
    > circuit.stim
stim analyze_errors \
    --decompose_errors \
    --fold_loops \
    --in circuit.stim \
    > error_model.dem
stim detect \
    --in circuit.stim \
    --shots $shots \
    --obs_out actual_obs_flips.01 \
    --obs_out_format 01 \
    --out detection_events.b8 \
    --out_format b8

for ((i=0; i<${#threads[@]}; i++ )); do
    thisThreads=${threads[$i]}
    mkdir -p log_${thisThreads}threads
    sbatch \
        --chdir=log_${thisThreads}threads \
        --output=../log_${thisThreads}threads.out \
        --cpus-per-task=$thisThreads \
        ~/PyMatchingSHMEM/scripts/profile_threads_call.sh $thisThreads $M
done
