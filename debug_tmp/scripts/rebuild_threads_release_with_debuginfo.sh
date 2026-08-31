#!/bin/bash
#SBATCH --job-name=rebuild_threads_release_g
#SBATCH --output=debug_tmp/out/rebuild_threads_release_g-%j.out
#SBATCH --error=debug_tmp/out/rebuild_threads_release_g-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB

# Adds -g to the Release build so VTune can attribute samples to source lines
# (not just function symbols) -- decode_shots._omp_fn.0 shows as a single flat
# self-time bucket without it because small helpers (wait_until_decode_done,
# wait_until_done, etc.) get fully inlined at -O3 with no debug info to unwind
# through. -g doesn't change codegen/-O3, only adds DWARF tables.

source ~/.bash_profile
conda activate pymatching
cd ~/PyMatchingSHMEM
cmake -DCMAKE_CXX_FLAGS_RELEASE="-O3 -g -DNDEBUG" -S . -B build_threads_release
cmake --build build_threads_release --target pymatching -j 16
