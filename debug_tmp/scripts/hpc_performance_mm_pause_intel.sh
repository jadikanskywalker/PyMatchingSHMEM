#!/bin/bash
#SBATCH --job-name=hpc_perf_mm_pause_intel
#SBATCH --output=debug_tmp/out/sbatch/hpc_perf_mm_pause_intel-%j.out
#SBATCH --error=debug_tmp/out/sbatch/hpc_perf_mm_pause_intel-%j.err
#SBATCH --partition=h100
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G

# hpc-performance collection for each thread count already covered by the hotspots sweep in
# run_d21_profile_mm_pause_intel -- vtune_profile.sh's topdown/hotspots report calls don't apply
# to hpc-performance's report schema (no Module column there), so this runs collect + summary
# report directly rather than going through vtune_profile.sh.

source ~/.bash_profile
conda activate pymatching

dir=~/PyMatchingSHMEM/run_d21_profile_mm_pause_intel
threads_build=~/PyMatchingSHMEM/build_threads_release/pymatching

export OMP_PLACES=cores
export OMP_PROC_BIND=close

for t in 16 32 64; do
    export OMP_NUM_THREADS=$t
    result_dir="${dir}/log_${t}threads/run/vtune/threads${t}_hpc"
    mkdir -p "$result_dir"

    vtune -collect hpc-performance \
          -knob enable-stack-collection=true \
          -r "${result_dir}/threads${t}_hpc" \
          --app-working-dir "${dir}/log_${t}threads" \
          -- $threads_build predict \
             --dem "${dir}/error_model.dem" \
             --in "${dir}/detection_events.b8" \
             --in_format b8 \
             --out "${dir}/log_${t}threads/predicted_obs_flips_hpc.01" \
             --out_format 01 \
             --rounds_per_partition 32 \
             --use_threads

    vtune -report summary -r "${result_dir}/threads${t}_hpc/threads${t}_hpc.vtune" \
          -report-output "${result_dir}/threads${t}_hpc_summary.txt"

    echo "Summary: ${result_dir}/threads${t}_hpc_summary.txt"
done
