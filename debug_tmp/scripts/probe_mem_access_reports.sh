#!/bin/bash
#SBATCH --job-name=probe_mem_access_reports
#SBATCH --output=debug_tmp/out/probe_mem_access_reports-%j.out
#SBATCH --error=debug_tmp/out/probe_mem_access_reports-%j.err
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
result=~/PyMatchingSHMEM/debug_tmp/out/probe_mem_access_result

export OMP_NUM_THREADS=16
export OMP_PLACES=cores
export OMP_PROC_BIND=close

vtune -collect memory-access \
      -knob analyze-mem-objects=true \
      -knob dram-bandwidth-limits=false \
      -r "$result" \
      -- $threads_build predict \
         --dem "${dir}/error_model.dem" \
         --in "${dir}/detection_events.b8" \
         --in_format b8 \
         --out ~/PyMatchingSHMEM/debug_tmp/out/probe_mem_access_out.01 \
         --out_format 01 \
         --rounds_per_partition 32 \
         --use_threads

echo "=== probing invalid report name for valid-list hint ==="
vtune -report bogus_report_name -r "$result" 2>&1

echo "=== summary report ==="
vtune -report summary -r "$result" -report-output=~/PyMatchingSHMEM/debug_tmp/out/probe_mem_summary.txt
echo "summary exit: $?"

echo "=== group-by=? for hotspots (memory objects?) ==="
vtune -report hotspots -r "$result" -group-by=? 2>&1 | head -60
