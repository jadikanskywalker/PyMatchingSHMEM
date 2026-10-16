#!/bin/bash
#SBATCH --job-name=test_pmu_zen4
#SBATCH --output=debug_tmp/out/sbatch/test_pmu_zen4-%j.out
#SBATCH --error=debug_tmp/out/sbatch/test_pmu_zen4-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

source ~/.bash_profile
conda activate pymatching

dir=~/PyMatchingSHMEM/run_d21_profile
threads_build=~/PyMatchingSHMEM/build_threads_release/pymatching

export OMP_NUM_THREADS=16
export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "=== memory-access test ==="
vtune -collect memory-access \
      -knob analyze-mem-objects=true \
      -r /tmp/test_mem_access \
      -- $threads_build predict \
         --dem "${dir}/error_model.dem" \
         --in "${dir}/detection_events.b8" \
         --in_format b8 \
         --out /tmp/test_mem_out.01 \
         --out_format 01 \
         --rounds_per_partition 32 \
         --use_threads
echo "memory-access exit: $?"

echo "=== uarch-exploration test ==="
vtune -collect uarch-exploration \
      -r /tmp/test_uarch \
      -- $threads_build predict \
         --dem "${dir}/error_model.dem" \
         --in "${dir}/detection_events.b8" \
         --in_format b8 \
         --out /tmp/test_uarch_out.01 \
         --out_format 01 \
         --rounds_per_partition 32 \
         --use_threads
echo "uarch-exploration exit: $?"
