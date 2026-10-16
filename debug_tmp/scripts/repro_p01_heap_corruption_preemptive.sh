#!/bin/bash
#SBATCH --job-name=repro_p01_heap_corruption_preemptive
#SBATCH --output=debug_tmp/out/sbatch/repro_p01_heap_corruption_preemptive-%j.out
#SBATCH --error=debug_tmp/out/sbatch/repro_p01_heap_corruption_preemptive-%j.err
#SBATCH --partition=zen4
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=64G

# Same repro as repro_p01_heap_corruption.sh, but with --extract_preemptively -- checking whether
# preemptive extraction changes the heap-corruption crash, and whether prune_stale_regions_matched_to_vb
# (mwpm.cc's shatter paths + decoding_unit.cc's post-hoc scan) ever fires on a plain ROUND-partitioning
# run with no SpecialTasks (LocalSeamTask/CrossRankTask) at all -- single-PE, no USE_SHMEM here, so any
# fire would have to come from the ordinary preemptive extraction-unit chain/back_divide_walk path, not
# from cross-rank or local-seam fusion.

source ~/.bash_profile
conda activate pymatching

threads_build=~/PyMatchingSHMEM/build_threads/pymatching
dem=~/PyMatchingSHMEM/debug_tmp/out/scorep_variance/p01/error_model.dem
det=~/PyMatchingSHMEM/debug_tmp/out/scorep_variance/p01/detection_events_1000.b8

outdir=~/PyMatchingSHMEM/debug_tmp/out/repro_p01_heap_corruption_preemptive
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
        --extract_preemptively \
        > run.out 2> run.err
    echo "${tag} exit: $?"
    cd ..
done
