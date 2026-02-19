#!/bin/bash

# conda activate pymatching
# module unload python
# module load python

# export SHMEM_OFI_PROVIDER=ofi_rxm

if [ $# -le 3 ]
  then
    echo "Args: [nthreads] [shots] [rounds] [M]"
    exit 1
else
    nthreads=$1
    shots=$2
    rounds=$(($3-1))
    M=$4
fi

ppn=2

if [ ! -d "run" ]
  then
    mkdir run
fi

cd run

rm *.out

hosts=$(srun hostname | sort | uniq | paste -sd, -)
# Function to create a hostfile with specified slots per host
create_hostfile() {
  local hostfile="hostfile.txt"
  # Clear the hostfile if it exists
  > "$hostfile"
  # Write each host and its slots to the hostfile
  for host in ${hosts//,/ }; do
    echo "$host slots=$ppn" >> "$hostfile"
  done
}
create_hostfile

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

# need to find d=21 lattice surgery circuit
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

# Run prediction
if [ $nthreads -gt 0 ]
  then
    export OMP_NUM_THREADS=$nthreads
    export OMP_PLACES=cores
    export OMP_PROC_BIND=true
fi

echo "Starting threads run..."
start_serial=$(date +%s)
~/PyMatchingSHMEM/build_threads_release/pymatching predict \
    --dem error_model.dem \
    --in detection_events.b8 \
    --in_format b8 \
    --out predicted_obs_flips__threads.01 \
    --out_format 01 \
    --rounds_per_partition $M \
    --use_threads \
    > log_threads.out
end_serial=$(date +%s)
serial_time=$((end_serial - start_serial))
echo "Threads run completed in $serial_time seconds."

echo Threads
echo correct predictions:
paste -d " " predicted_obs_flips__threads.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips__threads.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l
echo

# Parallel run with timing
# export ASAN_OPTIONS=detect_leaks=0
echo "Starting SHMEM run..."
start_parallel=$(date +%s)
oshrun  \
    --map-by ppr:${ppn}:node:pe=${OMP_NUM_THREADS} \
    --bind-to core \
    -hostfile hostfile.txt \
    ~/PyMatchingSHMEM/build_sos/pymatching predict \
    --dem error_model.dem \
    --in detection_events.b8 \
    --in_format b8 \
    --out predicted_obs_flips__shmem.01 \
    --out_format 01 \
    --rounds_per_partition $M \
    --use_threads \
    > log_shmem.out

end_parallel=$(date +%s)
parallel_time=$((end_parallel - start_parallel))
echo "SHMEM run completed in $parallel_time seconds."

# Check work
echo SHMEM
echo correct predictions:
paste -d " " predicted_obs_flips__shmem.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips__shmem.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l

echo
echo Shots with differring predictions:
awk 'NR==FNR{a[NR]=$0; n=NR; next} {
  if (FNR>n || $0!=a[FNR]) { print FNR-1; out=1 }
} END {
  if (n>FNR) { for (i=FNR+1;i<=n;i++) { print i-1; out=1 } }
  if (!out) print "no differences"
}' predicted_obs_flips__threads.01 predicted_obs_flips__shmem.01

rm hostfile.txt
cd ..