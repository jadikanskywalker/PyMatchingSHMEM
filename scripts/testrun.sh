#!/bin/bash

# conda activate pymatching
# module unload python
# module load python

# export SHMEM_OFI_PROVIDER=ofi_rxm

if [ $# -le 3 ]
  then
    echo "Args: [nthreads] [shots] [rounds] [M] ... [skip_stim]"
    exit 1
else
    nthreads=$1
    shots=$2
    rounds=$(($3-1))
    M=$4
fi

if [[ -n ${5+x} ]]; then
    build_circuit=false
else
    build_circuit=true
fi


serial_build=~/PyMatchingSHMEM/build/pymatching
threads_build=~/PyMatchingSHMEM/build_threads/pymatching
echo "serial_build:  $serial_build"
echo "threads_build: $threads_build"

if [ ! -d "run" ]
  then
    mkdir run
fi

cd run

if [ -d "out_parallel" ]
  then
    rm out_parallel -r
fi
if [ -d "out_serial" ]
  then
    rm out_serial -r
fi
if [ -d "out_frames" ]
  then
    rm out_frames -r
fi

if $build_circuit; then
    rm *.01 circuit.stim *.b8 *.dem
    python3 ../scripts/gen_multi_obs.py \
        --num_observables 4 \
        --rounds $rounds \
        --distance 7 \
        --after_clifford_depolarization 0.1 \
        --code repetition_code \
        --task memory \
        --num_surgery_gates 5 \
        --surgery_duration 5 \
        --circuit_out circuit.stim \
        > error_model.dem
    # Sample detection events FROM THE DEM (not the circuit) so seam errors fire
    stim sample_dem \
        --in error_model.dem \
        --shots $shots \
        --out detection_events.b8 \
        --out_format b8 \
        --obs_out actual_obs_flips.01 \
        --obs_out_format 01


    # stim gen \
    #     --rounds $rounds \
    #     --distance 5 \
    #     --after_clifford_depolarization 0.1 \
    #     --code repetition_code \
    #     --task memory \
    #     > circuit.stim
    # stim analyze_errors \
    #     --decompose_errors \
    #     --fold_loops \
    #     --in circuit.stim \
    #     > error_model.dem
    # stim detect \
    #     --in circuit.stim \
    #     --shots $shots \
    #     --obs_out actual_obs_flips.01 \
    #     --obs_out_format 01 \
    #     --out detection_events.b8 \
    #     --out_format b8

    # stim gen \
    #     --rounds=$rounds \
    #     --distance=7 \
    #     --after_clifford_depolarization=0.01 \
    #     --code surface_code \
    #     --task rotated_memory_x \
    #     > circuit.stim
    # stim analyze_errors \
    #     --decompose_errors \
    #     --fold_loops \
    #     --in circuit.stim \
    #     > error_model.dem
    # stim detect \
    #     --in circuit.stim \
    #     --shots $shots \
    #     --obs_out actual_obs_flips.01 \
    #     --obs_out_format 01 \
    #     --out detection_events.b8 \
    #     --out_format b8

    # echo "Starting serial run..."
    # start_serial=$(date +%s)
    # $serial_build predict \
    #     --dem error_model.dem \
    #     --in detection_events.b8 \
    #     --in_format b8 \
    #     --out predicted_obs_flips.01 \
    #     --out_format 01 \
    #     > log_serial.out
    # end_serial=$(date +%s)
    # serial_time=$((end_serial - start_serial))
    # echo "Serial run completed in $serial_time seconds."

    # echo Serial
    # echo correct predictions:
    # paste -d " " predicted_obs_flips.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
    # echo wrong predictions:
    # paste -d " " predicted_obs_flips.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l
    # echo
fi

# Run prediction
if [ $nthreads -gt 0 ]
  then
    export OMP_NUM_THREADS=$nthreads
    export OMP_PLACES=cores
    export OMP_PROC_BIND=true
fi

# Parallel run with timing
echo "Starting parallel run..."
start_parallel=$(date +%s)
$threads_build predict \
    --dem error_model.dem \
    --in detection_events.b8 \
    --in_format b8 \
    --out predicted_obs_flips__threads.01 \
    --out_format 01 \
    --rounds_per_partition $M \
    --seam_buffer_size $k \
    --task_division_strategy observable \
    --use_threads \
    --extraction_unit_size 4 \
    --extract_preemptively \
    --draw_frames \
    > log_threads.out 2>log_threads.err
end_parallel=$(date +%s)
parallel_time=$((end_parallel - start_parallel))
echo "Threads run completed in $parallel_time seconds."

# Check work
echo Threads
echo correct predictions:
paste -d " " predicted_obs_flips__threads.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips__threads.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l

# echo
# echo Shots with differring predictions:
# awk 'NR==FNR{a[NR]=$0; n=NR; next} {
#   if (FNR>n || $0!=a[FNR]) { print FNR-1; out=1 }
# } END {
#   if (n>FNR) { for (i=FNR+1;i<=n;i++) { print i-1; out=1 } }
#   if (!out) print "no differences"
# }' predicted_obs_flips.01 predicted_obs_flips__threads.01

rm hostfile.txt
cd ..