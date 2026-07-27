#!/bin/bash
#SBATCH --job-name=PyMatchingSHMEM
#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --partition=zen4
#SBATCH --time=03:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=1000GB

cd ~/PyMatchingSHMEM
source ~/.bashrc
conda activate pymatching
# module unload python
# module load python

# export FI_LOG_LEVEL=trace/
# export SHMEM_OFI_PROVIDER="verbs"
export FI_VERBS_DEVICE_NAME="mlx5_2"
# unset FI_VERBS_DEVICE_NAME
# export SHMEM_OFI_PROVIDER=shm

if [ $# -le 7 ]
  then
    echo "Args: [n] [ppn] [nthreads_shmem] [nthreads] [shots] [rounds] [M] [k]"
    exit 1
else
    n=$1
    ppn=$2
    nthreads_shmem=$3
    nthreads=$4
    shots=$5
    rounds=$(($6-1))
    M=$7
    k=$8
fi

if [ ! -d "run" ]
  then
    mkdir run
fi

cd run
rm *.out *.err

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
if [ -d "scorep_results" ]
  then
    rm -r scorep_results
fi

# # need to find d=21 lattice surgery circuit
# rm *.01 circuit.stim *.b8 *.dem
python3 ../scripts/gen_multi_obs.py \
    --num_observables 4 \
    --rounds $rounds \
    --distance 7 \
    --after_clifford_depolarization 0.1 \
    --code repetition_code \
    --task memory \
    --num_surgery_gates 3 \
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


# python scripts/gen_multi_obs.py \
#     --num_observables 24 \
#     --rounds 755 \
#     --distance 21 \
#     --after_clifford_depolarization 0.001 \
#     --code surface_code \
#     --task rotated_memory_x \
#     --surgery_preset 24obs \
#     --circuit_out circuit.stim \
#     > error_model.dem

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

# Run prediction
if [ $nthreads -gt 0 ]
  then
    export OMP_NUM_THREADS=$nthreads
    export OMP_PLACES=cores
    export OMP_PROC_BIND=true
fi

# echo "Starting serial run..."
# start_serial=$(date +%s)
# ~/PyMatching/build/pymatching predict \
#     --dem error_model.dem \
#     --in detection_events.b8 \
#     --in_format b8 \
#     --out predicted_obs_flips__threads.01 \
#     --out_format 01 \
#     > log_serial.out
#     # --rounds_per_partition $M \
#     # --use_threads \
#     # --draw_frames \
# end_serial=$(date +%s)
# serial_time=$((end_serial - start_serial))
# echo "Serial run completed in $serial_time seconds."

# echo Serial
# echo correct predictions:
# paste -d " " predicted_obs_flips__threads.01 actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
# echo wrong predictions:
# paste -d " " predicted_obs_flips__threads.01 actual_obs_flips.01 | grep "0 1\|1 0" | wc -l
# echo

# enable profiling
# export SCOREP_ENABLE_PROFILING=true
# export SCOREP_ENABLE_TRACING=false
# export SCOREP_EXPERIMENT_DIRECTORY=scorep_results
# export SCOREP_OVERWRITE_EXPERIMENT_DIRECTORY=true

# Parallel run with timing
echo "Starting SHMEM run..."
start_parallel=$(date +%s)
export OMP_NUM_THREADS=$nthreads_shmem
export SHMEM_SYMMETRIC_SIZE=16G
# export LIBFABRIC_DEBUG=0
$SWHOME/sos_1.5_scalable/bin/oshrun  \
    -n $n \
    --map-by ppr:$ppn:package:PE=$nthreads_shmem \
    --bind-to core \
    --report-bindings \
    ~/PyMatchingSHMEM/build_sos/pymatching predict \
    --dem error_model.dem \
    --in detection_events.b8 \
    --in_format b8 \
    --out predicted_obs_flips__shmem.01 \
    --out_format 01 \
    --rounds_per_partition $M \
    --obs_coors_included \
    --cross_rank_fusion_window_size $k \
    --task_division_strategy observable \
    --use_threads \
    --extraction_unit_size 2 \
    --draw_frames \
    > log_shmem.out 2>log_shmem.err


    # --dem ../testdems/error_model_48obs_d21_p001_640r.dem \
    # --in ../testdems/detection_events_48obs_d21_p001_640r_100s.b8 \
    # ~/PyMatchingSHMEM/scripts/scorep_wrapper.sh \
# gdb -q -batch -ex run -ex "thread apply all bt full" -ex "set logging file gdb.log" -ex 'info sharedlibrary libfabric' --args \
    # --dem ../testdems/error_model_9obs_d21_p001_672r.dem \
    # --in ../testdems/detection_events_9obs_d21_p001_672r_100s.b8 \

end_parallel=$(date +%s)
parallel_time=$((end_parallel - start_parallel))
echo "SHMEM run completed in $parallel_time seconds."

npes=$(($n * $ppn))

python3 ../scripts/combine_results.py predicted_obs_flips__shmem.01 $n

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
