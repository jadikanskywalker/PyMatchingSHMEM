#!/bin/bash

# conda activate pymatching
# module unload python
# module load python

# export SHMEM_OFI_PROVIDER=ofi_rxm

if [ $# -le 1 ]
  then
    echo "Args: [shots] [rounds]"
    exit 1
else
    shots=$1
    rounds=$(($2-1))
fi

if [ ! -d "run" ]
  then
    mkdir run
fi

cd run

rm *.01 circuit.stim *.b8 *.dem

if [ -d "out_serial" ]
  then
    rm out_serial -r
fi
if [ -d "out_frames" ]
  then
    rm out_frames -r
fi

stim gen \
    --rounds=$rounds \
    --distance=21 \
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

echo "Starting serial run..."
start_serial=$(date +%s)
~/PyMatchingSHMEM/build/pymatching predict \
    --dem error_model.dem \
    --in detection_events.b8 \
    --in_format b8 \
    --out predicted_obs_flips__without_shmem.01 \
    --out_format 01 \
    > log_serial.out
end_serial=$(date +%s)
serial_time=$((end_serial - start_serial))
echo "Serial run completed in $serial_time seconds."

echo Serial
echo correct predictions:
paste -d " " predicted_obs_flips__without_shmem.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips__without_shmem.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l
echo

cd ..