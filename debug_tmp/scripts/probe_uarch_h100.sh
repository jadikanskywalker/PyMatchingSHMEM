#!/bin/bash
#SBATCH --job-name=probe_uarch_h100
#SBATCH --output=debug_tmp/out/probe_uarch_h100-%j.out
#SBATCH --error=debug_tmp/out/probe_uarch_h100-%j.err
#SBATCH --partition=h100
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

source ~/.bash_profile
conda activate pymatching

dir=~/PyMatchingSHMEM/run_d21_profile
threads_build=~/PyMatchingSHMEM/build_threads_release/pymatching
result=~/PyMatchingSHMEM/debug_tmp/out/probe_uarch_result

export OMP_NUM_THREADS=16
export OMP_PLACES=cores
export OMP_PROC_BIND=close

vtune -collect uarch-exploration \
      -r "$result" \
      -- $threads_build predict \
         --dem "${dir}/error_model.dem" \
         --in "${dir}/detection_events.b8" \
         --in_format b8 \
         --out ~/PyMatchingSHMEM/debug_tmp/out/probe_uarch_out.01 \
         --out_format 01 \
         --rounds_per_partition 32 \
         --use_threads
echo "collect exit: $?"

vtune -report summary -r "$result" -report-output=~/PyMatchingSHMEM/debug_tmp/out/probe_uarch_summary.txt
echo "summary exit: $?"

echo "=== group-by=? for hotspots ==="
vtune -report hotspots -r "$result" -group-by=? 2>&1
