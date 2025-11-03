#!/bin/bash

# module unload python
# module load python

# module load intel-oneapi-compilers
# module load intel-oneapi-mpi

rm -rf logs/vtune_result.*
rm logs/vtune_test_topdown.csv

export SHMEM_OFI_PROVIDER=ofi_rxm

if [ $# -eq 0 ]
  then
    ppn=1
    M=10
elif [ $# -eq 1 ]
  then
    ppn=$1
    M=10
else
    ppn=$1
    M=$2
fi

rm *.01 circuit.stim *.b8 *.dem

stim gen \
    --rounds=100 \
    --distance=5 \
    --after_clifford_depolarization=0.003 \
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
    --shots 1000 \
    --obs_out actual_obs_flips.01 \
    --obs_out_format 01 \
    --out detection_events.b8 \
    --out_format b8

# Run prediction
oshrun --hostfile hostfile.txt -N $ppn \
  vtune -collect hotspots -knob enable-stack-collection=true -r logs/vtune_result \
  ./pymatching predict \
    --dem error_model.dem \
    --in detection_events.b8 \
    --in_format b8 \
    --out predicted_obs_flips.01 \
    --out_format 01 \
    --rounds_per_partition $M

vtune -report top-down -r logs/vtune_result.rpg-* \
      -call-stack-mode all -column="CPU Time:Self","Module" -filter "Function Stack"  \
      --format csv -csv-delimiter comma -report-output logs/vtune_report_topdown.csv
      

# Check work
echo correct predictions:
paste -d " " predicted_obs_flips.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l
