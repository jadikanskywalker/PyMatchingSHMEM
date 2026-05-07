#!/bin/bash
#SBATCH --job-name=bench_threads
#SBATCH --output=bench_threads-%j.out
#SBATCH --partition=zen4
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --sockets=1
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=8G

if [ $# -le 3 ]
  then
    echo "Args: [d] [p_dec] [shots] [dir] ... [skip_stim]"
    exit 1
else
    d=$1
    p_dec=$2
    shots=$3
    dir=$4
fi

if [[ -n ${5+x} ]]; then
    build_circuit=false
else
    build_circuit=true
fi

source ~/.bash_profile
conda activate pymatching

if [ ! -d $dir ]
  then
    mkdir $dir
fi

cd $dir
pwd > bench.out

rounds=8192
threads=(2 4 8 16 32)
# threads=(32)
# M=(64 32 16 8)
M=(8 16 32 64)
# M=(128)
# threads=(16)
# M=(8)
code=surface_code
task=rotated_memory_x
p=0.$p_dec

serial_build=~/PyMatchingSHMEM/build/pymatching
threads_build=~/PyMatchingSHMEM/build_threads_release/pymatching
echo "serial_build:  $serial_build" >> bench.out
echo "threads_build: $threads_build"

if $build_circuit; then
    stim gen \
        --rounds=$(($rounds-1)) \
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

    # rerun serial for new circuit
    # start_serial=$(date +%s)
    # $serial_build predict \
    #     --dem error_model.dem \
    #     --in detection_events.b8 \
    #     --in_format b8 \
    #     --out predicted_obs_flips.01 \
    #     --out_format 01 \
    #     > log_0.out
    # end_serial=$(date +%s)
    # serial_time=$((end_serial - start_serial))
    # echo "0: $serial_time seconds"
fi

for ((m=0; m<${#M[@]}; m++ )); do
    thisM=${M[$m]}
    echo "----------"
    echo "M: $thisM"
    for ((i=0; i<${#threads[@]}; i++ )); do
        thisThreads=${threads[$i]}

        sbatch \
            --output=log_M${thisM}_${thisThreads}threads.out \
            ~/PyMatchingSHMEM/scripts/benchmark_call.sh $thisM $thisThreads
    done
done

start_serial=$(date +%s)
$serial_build predict \
    --dem error_model.dem \
    --in detection_events.b8 \
    --in_format b8 \
    --out predicted_obs_flips.01 \
    --out_format 01 \
    --num_repeats 10 \
    > log_0.out
end_serial=$(date +%s)
serial_time=$((end_serial - start_serial))
echo "0: $serial_time seconds"

