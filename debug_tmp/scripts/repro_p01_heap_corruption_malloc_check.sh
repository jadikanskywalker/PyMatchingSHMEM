#!/bin/bash
#SBATCH --job-name=repro_p01_heap_corruption_malloc_check
#SBATCH --output=debug_tmp/out/sbatch/repro_p01_heap_corruption_malloc_check-%j.out
#SBATCH --error=debug_tmp/out/sbatch/repro_p01_heap_corruption_malloc_check-%j.err
#SBATCH --partition=zen4
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=64G

# Same repro as repro_p01_heap_corruption.sh, plus glibc's own MALLOC_CHECK_=3 alongside ASan.
# MALLOC_CHECK_=3 makes glibc abort with its own diagnostic the instant it detects corrupted heap
# metadata (double-free, invalid free, heap overflow) in its own allocator bookkeeping.

source ~/.bash_profile
conda activate pymatching

export MALLOC_CHECK_=3

threads_build=~/PyMatchingSHMEM/build_threads/pymatching
dem=~/PyMatchingSHMEM/debug_tmp/out/scorep_variance/p01/error_model.dem
det=~/PyMatchingSHMEM/debug_tmp/out/scorep_variance/p01/detection_events_1000.b8

outdir=~/PyMatchingSHMEM/debug_tmp/out/repro_p01_heap_corruption_malloc_check
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
