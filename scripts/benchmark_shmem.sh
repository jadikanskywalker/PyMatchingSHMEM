#!/bin/bash
#SBATCH --job-name=pymatching
#SBATCH --output=log001.out
#SBATCH --error=log001.err
#SBATCH --partition=zen4
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --exclusive

if [ $# -le 2 ]
  then
    echo "Args: [d] [p_dec] [shots] ... [skip_stim]"
    exit 1
else
    d=$1
    p_dec=$2
    shots=$3
fi

if [[ -n ${4+x} ]]; then
    build_circuit=false
else
    build_circuit=true
fi

source ~/.bash_profile
conda activate pymatching

if [ ! -d "bench$p_dec" ]
  then
    mkdir bench$p_dec
fi

cd bench$p_dec

rounds=32768
threads=(8 16 32 64 128)
M=(32 64)
# shmem_n=(2 4)
shmem_n2_ppn=(2 2 2 2 2 1)
shmem_n2_threads=(4 8 16 32 64 128)
shmem_n4_ppn=(4 4 4 4 2 2)
shmem_n4_threads=(4 8 16 32 32 64)
code=surface_code
task=rotated_memory_x
p=0.$p_dec

serial_build=~/PyMatchingSHMEM/build/pymatching
threads_build=~/PyMatchingSHMEM/build_threads_release/pymatching
shmem_build=~/PyMatchingSHMEM/build_sos/pymatching
echo "serial_build:  $serial_build" >> bench.out
echo "threads_build: $threads_build" >> bench.out
echo "shmem_build: $shmem_build" >> bench.out

# if $build_circuit; then
#     stim gen \
#         --rounds=$(($rounds-1)) \
#         --distance=$d \
#         --after_clifford_depolarization=$p \
#         --code $code \
#         --task $task \
#         > circuit.stim
#     stim analyze_errors \
#         --decompose_errors \
#         --fold_loops \
#         --in circuit.stim \
#         > error_model.dem
#     stim detect \
#         --in circuit.stim \
#         --shots $shots \
#         --obs_out actual_obs_flips.01 \
#         --obs_out_format 01 \
#         --out detection_events.b8 \
#         --out_format b8
# fi

start_serial=$(date +%s)
$serial_build predict \
    --dem error_model.dem \
    --in detection_events.b8 \
    --in_format b8 \
    --out predicted_obs_flips.01 \
    --out_format 01 \
    > log_0.out
end_serial=$(date +%s)
serial_time=$((end_serial - start_serial))
echo "----------" >> bench.out
echo "0: $serial_time seconds" >> bench.out

export OMP_PLACES=cores
export OMP_PROC_BIND=true

export SHMEM_SYMMETRIC_SIZE=4G

for ((m=0; m<${#M[@]}; m++ )); do
    thisM=${M[$m]}
    echo "----------" >> bench.out
    echo "M: $thisM" >> bench.out
    echo "  THREADS:" >> bench.out
    for ((i=0; i<${#threads[@]}; i++ )); do
        thisThreads=${threads[$i]}

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
            > log_M${thisM}_${thisThreads}threads.out
        end_parallel=$(date +%s)
        parallel_time=$((end_parallel - start_parallel))
        echo "    Threads=$thisThreads: $parallel_time seconds" >> bench.out
    done
    # echo "  SHMEM:" >> bench.out
    # for ((i=0; i<${#shmem_n2_ppn[@]}; i++ )); do
    #     thisPPN=${shmem_n2_ppn[$i]}
    #     thisThreads=${shmem_n2_threads[$i]}

    #     export OMP_NUM_THREADS=$thisThreads

    #     start_parallel=$(date +%s)
    #     oshrun  \
    #         -n 2 \
    #         --map-by ppr:$thisPPN:node:PE=$thisThreads \
    #         --bind-to core \
    #         --report-bindings \
    #         $shmem_build predict \
    #             --dem error_model.dem \
    #             --in detection_events.b8 \
    #             --in_format b8 \
    #             --out predicted_obs_flips.01 \
    #             --out_format 01 \
    #             --rounds_per_partition $thisM \
    #             --cross_rank_fusion_window_size 1 \
    #             --use_threads \
    #             > log_M${thisM}_n2_ppn${thisPPN}_${thisThreads}threads_k1.out
    #     end_parallel=$(date +%s)
    #     parallel_time=$((end_parallel - start_parallel))
    #     echo "  N=2 PPN=$thisPPN Threads=$thisThreads: $parallel_time seconds" >> bench.out
    # done
    # for ((i=0; i<${#shmem_n4_ppn[@]}; i++ )); do
    #     thisPPN=${shmem_n4_ppn[$i]}
    #     thisThreads=${shmem_n4_threads[$i]}

    #     export OMP_NUM_THREADS=$thisThreads

    #     start_parallel=$(date +%s)
    #     oshrun  \
    #         -n 4 \
    #         --map-by ppr:$thisPPN:node:PE=$thisThreads \
    #         --bind-to core \
    #         --report-bindings \
    #         $shmem_build predict \
    #             --dem error_model.dem \
    #             --in detection_events.b8 \
    #             --in_format b8 \
    #             --out predicted_obs_flips.01 \
    #             --out_format 01 \
    #             --rounds_per_partition $thisM \
    #             --cross_rank_fusion_window_size 1 \
    #             --use_threads \
    #             > log_M${thisM}_n4_ppn${thisPPN}_${thisThreads}threads_k1.out
    #     end_parallel=$(date +%s)
    #     parallel_time=$((end_parallel - start_parallel))
    #     echo "  N=4 PPN=$thisPPN Threads=$thisThreads: $parallel_time seconds" >> bench.out
    # done
done