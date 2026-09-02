#!/bin/bash
#SBATCH --job-name=hpc_perf_preemptive_L16_intel_64t
#SBATCH --output=debug_tmp/out/hpc_perf_preemptive_L16_intel_64t-%j.out
#SBATCH --error=debug_tmp/out/hpc_perf_preemptive_L16_intel_64t-%j.err
#SBATCH --partition=h100
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G

# hpc-performance at 64 threads for preemptive extraction + extraction_unit_size 16, since
# benchmarking showed L=16 helps more at higher thread counts and preemptive gives a steady
# marginal win over non-preemptive at those L values.

source ~/.bash_profile
conda activate pymatching

dir=~/PyMatchingSHMEM/run_d21_profile_extraction_queue_preemptive_L16_intel_64t
threads_build=~/PyMatchingSHMEM/build_threads_release/pymatching
mkdir -p "$dir"
cp ~/PyMatchingSHMEM/run_d21_profile/circuit.stim ~/PyMatchingSHMEM/run_d21_profile/error_model.dem \
   ~/PyMatchingSHMEM/run_d21_profile/detection_events.b8 ~/PyMatchingSHMEM/run_d21_profile/actual_obs_flips.01 \
   "$dir/"

export OMP_NUM_THREADS=64
export OMP_PLACES=cores
export OMP_PROC_BIND=close

result_dir="${dir}/log_64threads/run/vtune/threads64_hpc"
mkdir -p "$result_dir"

vtune -collect hpc-performance \
      -knob enable-stack-collection=true \
      -r "${result_dir}/threads64_hpc" \
      --app-working-dir "${dir}/log_64threads" \
      -- $threads_build predict \
         --dem "${dir}/error_model.dem" \
         --in "${dir}/detection_events.b8" \
         --in_format b8 \
         --out "${dir}/log_64threads/predicted_obs_flips_hpc.01" \
         --out_format 01 \
         --rounds_per_partition 32 \
         --use_threads \
         --extract_preemptively \
         --extraction_unit_size 16

vtune -report summary -r "${result_dir}/threads64_hpc/threads64_hpc.vtune" \
      -report-output "${result_dir}/threads64_hpc_summary.txt"

echo "correct predictions:"
paste -d " " "${dir}/log_64threads/predicted_obs_flips_hpc.01" "${dir}/actual_obs_flips.01" | grep "1 1\|0 0" | wc -l
echo "wrong predictions:"
paste -d " " "${dir}/log_64threads/predicted_obs_flips_hpc.01" "${dir}/actual_obs_flips.01" | grep "0 1\|1 0" | wc -l

echo "Summary: ${result_dir}/threads64_hpc_summary.txt"
