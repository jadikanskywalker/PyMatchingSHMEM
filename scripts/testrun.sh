#!/bin/bash

# conda activate pymatching
# module unload python
# module load python

# export SHMEM_OFI_PROVIDER=ofi_rxm

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

stim gen \
    --rounds=$rounds \
    --distance=20 \
    --after_clifford_depolarization=0.03 \
    --code repetition_code \
    --task memory \
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

if [ $ppn -le 0 ]
  then
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

echo Serial
echo correct predictions:
paste -d " " predicted_obs_flips__without_shmem.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips__without_shmem.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l
echo

# Run prediction
if [ $nthreads -gt 0 ]
  then
    export OMP_NUM_THREADS=$nthreads
fi

if [ $ppn -le 0 ]
  then
    ~/PyMatchingSHMEM/build_threads/pymatching predict \
        --dem error_model.dem \
        --in detection_events.b8 \
        --in_format b8 \
        --out predicted_obs_flips__with_shmem.01 \
        --out_format 01 \
        --rounds_per_partition $M \
        --parallel \
        --draw_frames \
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