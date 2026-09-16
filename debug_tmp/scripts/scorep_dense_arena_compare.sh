#!/bin/bash
#SBATCH --job-name=scorep_dense_arena_compare
#SBATCH --output=debug_tmp/out/sbatch/scorep_dense_arena_compare-%j.out
#SBATCH --error=debug_tmp/out/sbatch/scorep_dense_arena_compare-%j.err
#SBATCH --partition=h100
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=64G

# Step B of the DenseArena-port plan (plans/snappy-toasting-bee.md): the actual cache-locality
# hypothesis test. Compares "Leaf Decode"'s Score-P cache-miss rate (cache-misses/cache-references)
# between the baseline Arena build (build_threads_scorep) and the new DenseArena build
# (build_threads_scorep_dense, ENABLE_DENSE_REGION_ARENA=ON) at matched thread counts. Step A
# already confirmed both builds produce byte-identical decode predictions, so any delta here is
# attributable purely to GraphFillRegion's memory layout.
#
# Reuses the same run_d21_profile_extraction_queue_intel workload (d=21, p=0.001, rounds=10000,
# M=32) as this session's earlier scorep_memory_metrics.sh, for direct comparability with the
# already-analyzed 86-89% Leaf Decode cache-miss-rate baseline on this same h100/Intel hardware.

source ~/.bash_profile
conda activate pymatching

dir=~/PyMatchingSHMEM/run_d21_profile_extraction_queue_intel
arena_build=~/PyMatchingSHMEM/build_threads_scorep/pymatching
dense_build=~/PyMatchingSHMEM/build_threads_scorep_dense/pymatching
outdir=~/PyMatchingSHMEM/debug_tmp/out/scorep_dense_arena_compare
mkdir -p "$outdir"
cd "$outdir"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export SCOREP_METRIC_PERF="cache-misses,cache-references"
export SCOREP_ENABLE_PROFILING=true
export SCOREP_ENABLE_TRACING=false
export SCOREP_TOTAL_MEMORY=1G

for t in 16 32 64; do
    export OMP_NUM_THREADS=$t
    for variant in arena dense; do
        if [ "$variant" == "arena" ]; then
            build=$arena_build
        else
            build=$dense_build
        fi
        tag="${variant}_threads${t}"
        mkdir -p "$tag"
        cd "$tag"
        export SCOREP_EXPERIMENT_DIRECTORY=scorep_results

        "$build" predict \
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
done
