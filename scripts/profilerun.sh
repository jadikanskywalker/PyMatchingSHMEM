#!/bin/bash

# module unload python
# module load python

# module load intel-oneapi-compilers
# module load intel-oneapi-mpi

rm -rf logs/vtune_result.*
rm logs/vtune_test_topdown.csv

export SHMEM_OFI_PROVIDER=ofi_rxm

if [ $# -le 4 ]
  then
    echo "Args: [ppn] [nthreads] [shots] [rounds] [M]"
    exit 1
else
    ppn=$1
    nthreads=$2
    shots=$3
    rounds=$4
    M=$5
fi

if [ ! -d "run" ]
  then
    mkdir run
fi

cd run

rm *.01 circuit.stim *.b8 *.dem

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


stim gen \
    --rounds=$rounds \
    --distance=7 \
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

# Run prediction
# oshrun --hostfile hostfile.txt -N $ppn \
#   vtune -collect hotspots -knob enable-stack-collection=true -r logs/vtune_result \
#   ./pymatching predict \
#     --dem error_model.dem \
#     --in detection_events.b8 \
#     --in_format b8 \
#     --out predicted_obs_flips.01 \
#     --out_format 01 \
#     --rounds_per_partition $M


echo "Starting serial run..."
start_serial=$(date +%s)
if [ $ppn -le 0 ]
  then
    vtune -collect hotspots -knob enable-stack-collection=true -r logs/serial/vtune_result \
      ~/PyMatchingSHMEM/build_threads/pymatching predict \
        --dem error_model.dem \
        --in detection_events.b8 \
        --in_format b8 \
        --out predicted_obs_flips__without_shmem.01 \
        --out_format 01 \
        --rounds_per_partition $M \
        > log_serial.out
else
    oshrun --hostfile hostfile.txt -N $ppn \
        ~/PyMatchingSHMEM/build_osss/pymatching predict \
        --dem error_model.dem \
        --in detection_events.b8 \
        --in_format b8 \
        --out predicted_obs_flips__without_shmem.01 \
        --out_format 01 \
        --rounds_per_partition $M > log_serial.out
fi
end_serial=$(date +%s)
serial_time=$((end_serial - start_serial))
echo "Serial run completed in $serial_time seconds."

echo Serial
echo correct predictions:
paste -d " " predicted_obs_flips__without_shmem.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips__without_shmem.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l
echo

vtune -report top-down -r logs/serial/vtune_result/vtune_result.vtune \
      -call-stack-mode all -column="CPU Time:Self","Module" -filter "Function Stack"  \
      --format csv -csv-delimiter comma -report-output logs/serial/vtune_report_topdown.csv
# mat b8

# Run prediction
if [ $nthreads -gt 0 ]
  then
    export OMP_NUM_THREADS=$nthreads
fi

# Parallel run with timing
echo "Starting parallel run..."
start_parallel=$(date +%s)
if [ $ppn -le 0 ]
  then
    vtune -collect hotspots -knob enable-stack-collection=true -r logs/parallel/vtune_result \
      ~/PyMatchingSHMEM/build_threads/pymatching predict \
        --dem error_model.dem \
        --in detection_events.b8 \
        --in_format b8 \
        --out predicted_obs_flips__with_shmem.01 \
        --out_format 01 \
        --rounds_per_partition $M \
        --parallel \
        > log_parallel.out
else
    oshrun --hostfile hostfile.txt -N $ppn \
      ~/PyMatchingSHMEM/build_osss/pymatching predict \
        --dem error_model.dem \
        --in detection_events.b8 \
        --in_format b8 \
        --out predicted_obs_flips__with_shmem.01 \
        --out_format 01 \
        --rounds_per_partition $M \
        --parallel \
        > log_parallel.out
fi
end_parallel=$(date +%s)
parallel_time=$((end_parallel - start_parallel))
echo "Parallel run completed in $parallel_time seconds."

vtune -report top-down -r logs/parallel/vtune_result/vtune_result.vtune \
      -call-stack-mode all -column="CPU Time:Self","Module" -filter "Function Stack"  \
      --format csv -csv-delimiter comma -report-output logs/parallel/vtune_report_topdown.csv
# mat b8

# Check work
echo Parallel
echo correct predictions:
paste -d " " predicted_obs_flips__with_shmem.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips__with_shmem.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l

echo
echo Shots with differring predictions:
awk 'NR==FNR{a[NR]=$0; n=NR; next} {
  if (FNR>n || $0!=a[FNR]) { print FNR-1; out=1 }
} END {
  if (n>FNR) { for (i=FNR+1;i<=n;i++) { print i-1; out=1 } }
  if (!out) print "no differences"
}' predicted_obs_flips__without_shmem.01 predicted_obs_flips__with_shmem.01

rm hostfile.txt
cd ..