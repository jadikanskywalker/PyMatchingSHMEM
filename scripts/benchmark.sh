#!/bin/bash
#SBATCH --job-name=pymatching001
#SBATCH --output=log001.out
#SBATCH --error=log001.err
#SBATCH --partition=h100
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G

if [ $# -le 2 ]
  then
    echo "Args: [d] [p_dec] [shots]"
    exit 1
else
    d=$1
    p_dec=$2
    shots=$3
fi

source ~/.bash_profile
conda activate pymatching

if [ ! -d "bench$p_dec" ]
  then
    mkdir bench$p_dec
fi

cd bench$p_dec

rm *.01 circuit.stim *.b8 *.dem *.out

rounds=128
threads=(2 4 8 16)
M=(64 32 16 8)
code=surface_code
task=rotated_memory_x
p=0.$p_dec

serial_build=~/PyMatchingSHMEM/build/pymatching
threads_build=~/PyMatchingSHMEM/build_threads/pymatching

stim gen \
    --rounds=$rounds \
    --distance=$d \
    --after_clifford_depolarization=$p \
    --code $code \
    --task $task \
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

# export OMP_NUM_THREADS=$i

start_parallel=$(date +%s)
$serial_build predict \
    --dem error_model.dem \
    --in detection_events.b8 \
    --in_format b8 \
    --out predicted_obs_flips.01 \
    --out_format 01 \
    > log_0.out
end_parallel=$(date +%s)
parallel_time=$((end_parallel - start_parallel))
echo "0: $parallel_time seconds"

export OMP_PLACES=cores
export OMP_PROC_BIND=close

for ((i=0; i<${#threads[@]}; i++ )); do
    # stim detect \
    #     --in circuit.stim \
    #     --shots $shots \
    #     --obs_out actual_obs_flips.01 \
    #     --obs_out_format 01 \
    #     --out detection_events.b8 \
    #     --out_format b8

    thisThreads=${threads[$i]}
    thisM=${M[$i]}

    export OMP_NUM_THREADS=$thisThreads

    start_parallel=$(date +%s)
    $threads_build predict \
        --dem error_model.dem \
        --in detection_events.b8 \
        --in_format b8 \
        --out predicted_obs_flips.01 \
        --out_format 01 \
        --rounds_per_partition $thisM \
        --use_threads \
        > log_$thisThreads.out
    end_parallel=$(date +%s)
    parallel_time=$((end_parallel - start_parallel))
    echo "$thisThreads: $parallel_time seconds"
done