#!/bin/bash
#SBATCH --job-name=profile_obs
#SBATCH --output=profile-%j.out
#SBATCH --partition=zen4
#SBATCH --time=05:00:00
#SBATCH --mem=1000GB

echo $SLURM_JOB_ID

cd ~/PyMatchingSHMEM
source ~/.bash_profile
conda activate pymatching

# export SHMEM_OFI_PROVIDER="verbs"
export FI_VERBS_DEVICE_NAME="mlx5_2"

if [ $# -le 9 ]
  then
    echo "Args: [n] [pps] [package/node] [nthreads_shmem] [M] [k] [dir] [dem] [det] [flips]"
    exit 1
else
    n=$1
    pps=$2
    perwhat=$3
    nthreads_shmem=$4
    M=$5
    k=$6
    dir=$7
    dem=$8
    det=$9
    flips=${10}
fi

if [ ! -d $dir ]
  then
    mkdir $dir
fi

cd $dir

rm scorep_results* -r
rm out_parallel -r
rm *.out *.err *.01 *.txt

# # need to find d=21 lattice surgery circuit
# rm *.01 circuit.stim *.b8 *.dem
# python3 ../scripts/gen_multi_obs.py \
#     --num_observables 2 \
#     --rounds 230 \
#     --distance 21 \
#     --after_clifford_depolarization 0.001 \
#     --code surface_code \
#     --task rotated_memory_x \
#     --surgery_spec "0,1,0,21;0,1,63,21;0,1,126,21;0,1,189,21"  \
#     --circuit_out circuit.stim \
#     > error_model_4obs_d21_p001_231r.dem
# # Sample detection events FROM THE DEM (not the circuit) so seam errors fire
# stim sample_dem \
#     --in error_model_4obs_d21_p001_231r.dem \
#     --shots 250 \
#     --out detection_events_4obs_d21_p001_231r_250s.b8 \
#     --out_format b8 \
#     --obs_out actual_obs_flips_4obs_d21_p001_231r_250s.01 \
#     --obs_out_format 01

# Run prediction
export OMP_NUM_THREADS=$nthreads
export OMP_PLACES=cores
export OMP_PROC_BIND=true


# enable profiling
export SCOREP_ENABLE_PROFILING=true
export SCOREP_ENABLE_TRACING=false
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
    --map-by ppr:$pps:$perwhat:PE=$nthreads_shmem \
    --bind-to core \
    --report-bindings \
    ~/PyMatchingSHMEM/scripts/scorep_wrapper.sh \
    ~/PyMatchingSHMEM/build_sos_profile/pymatching predict \
    --dem ../testdems/$dem \
    --in ../testdems/$det \
    --in_format b8 \
    --out predicted_obs_flips__shmem.01 \
    --out_format 01 \
    --rounds_per_partition $M \
    --obs_coors_included \
    --cross_rank_fusion_window_size $k \
    --task_division_strategy observable \
    --use_threads \
    --num_repeats 10 \
    &> log_shmem.out

end_parallel=$(date +%s)
parallel_time=$((end_parallel - start_parallel))
echo "SHMEM run completed in $parallel_time seconds."

python3 ../scripts/combine_results.py predicted_obs_flips__shmem.01 $n

# Check work
echo SHMEM
echo correct predictions:
paste -d " " predicted_obs_flips__shmem.01 ../testdems/$flips | grep "1 1\|0 0" | wc -l
echo wrong predictions:
paste -d " " predicted_obs_flips__shmem.01 ../testdems/$flips | grep "0 1\|1 0" | wc -l

cd ..
