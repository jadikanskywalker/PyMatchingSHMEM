#!/bin/bash
#SBATCH --job-name=scorep_memory_metrics
#SBATCH --output=debug_tmp/out/sbatch/scorep_memory_metrics-%j.out
#SBATCH --error=debug_tmp/out/sbatch/scorep_memory_metrics-%j.err
#SBATCH --partition=h100
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=64G

# Per-region hardware-counter breakdown via Score-P's native SCOREP_METRIC_PERF (no PAPI, no
# VTune driver needed -- confirmed available on this cluster's Score-P 9.4). Attaches cache-misses/
# cache-references to local_decoding/leaf_decode/solution_extraction/extraction_slot_wait, to
# directly test whether DRAM-bound stalls (37.8% Memory Bound / 11.1% Cache Bound / 29.9% DRAM
# Bound in the existing run_d21_profile_extraction_queue_intel hpc-performance data) concentrate in
# leaf decode (blossom-matching itself) vs. extraction vs. the wait loop's own atomic traffic.
#
# Reuses the exact existing error_model.dem/detection_events.b8 from
# run_d21_profile_extraction_queue_intel (d=21, p=0.001, rounds=10000, M=32) rather than
# regenerating, so results are directly comparable to that prior hpc-performance sweep.
#
# --exclusive rather than --sockets=1: this is specifically the memory-bandwidth-sensitive job of
# the three, so ruling out any other job's LLC/DRAM contention on the same node matters most here.

source ~/.bash_profile
conda activate pymatching

dir=~/PyMatchingSHMEM/run_d21_profile_extraction_queue_intel
scorep_build=~/PyMatchingSHMEM/build_threads_scorep/pymatching
outdir=~/PyMatchingSHMEM/debug_tmp/out/scorep_memory_metrics
mkdir -p $outdir
cd $outdir

export OMP_PLACES=cores
export OMP_PROC_BIND=close
# Generic, cross-platform perf event aliases (not vendor-specific uarch event names) -- this job
# was written targeting h100's Intel Sapphire Rapids nodes for direct comparability with the
# existing Intel hpc-performance data; if per-region breakdown needs vendor-specific events later,
# confirm exact names via `perf list` on the *target* node first (AMD/Intel PMU event names differ).
export SCOREP_METRIC_PERF="cache-misses,cache-references"
export SCOREP_ENABLE_PROFILING=true
export SCOREP_ENABLE_TRACING=false

for t in 16 32 64; do
    export OMP_NUM_THREADS=$t
    tag=threads${t}
    mkdir -p $tag
    cd $tag
    export SCOREP_EXPERIMENT_DIRECTORY=scorep_results

    $scorep_build predict \
        --dem "${dir}/error_model.dem" \
        --in "${dir}/detection_events.b8" \
        --in_format b8 \
        --out predicted_obs_flips.01 \
        --out_format 01 \
        --rounds_per_partition 32 \
        --use_threads \
        > run.out 2> run.err
    echo "${tag} exit: $?"
    cd ..
done
