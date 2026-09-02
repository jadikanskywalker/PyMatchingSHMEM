#!/bin/bash
#SBATCH --job-name=repro_p01_heap_corruption
#SBATCH --output=debug_tmp/out/repro_p01_heap_corruption-%j.out
#SBATCH --error=debug_tmp/out/repro_p01_heap_corruption-%j.err
#SBATCH --partition=zen4
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=64G

# Repro attempt for the glibc heap-corruption abort seen in the p=0.01 Score-P profile-mode sweep
# (16/32/64 threads all aborted with different malloc-corruption signatures; 8 threads succeeded).
# Uses the plain ASan-instrumented build_threads (no Score-P instrumentation at all), against the
# exact same p=0.01 DEM/detection-events already generated for that sweep, to determine whether
# this is a real pre-existing bug in the decoder or something specific to the Score-P build.

source ~/.bash_profile
conda activate pymatching

threads_build=~/PyMatchingSHMEM/build_threads/pymatching
dem=~/PyMatchingSHMEM/debug_tmp/out/scorep_variance/p01/error_model.dem
det=~/PyMatchingSHMEM/debug_tmp/out/scorep_variance/p01/detection_events_1000.b8

outdir=~/PyMatchingSHMEM/debug_tmp/out/repro_p01_heap_corruption
mkdir -p $outdir
cd $outdir

for t in 16 32 64; do
    export OMP_NUM_THREADS=$t
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close
    tag=threads${t}
    mkdir -p $tag
    cd $tag
    $threads_build predict \
        --dem "$dem" \
        --in "$det" \
        --in_format b8 \
        --out predicted_obs_flips.01 \
        --out_format 01 \
        --rounds_per_partition 32 \
        --use_threads \
        > run.out 2> run.err
    echo "${tag} exit: $?"
    cd ..
done
