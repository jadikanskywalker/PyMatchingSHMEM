#!/bin/bash
#SBATCH --job-name=rebuild_serial_pymatching
#SBATCH --output=debug_tmp/out/rebuild_serial_pymatching-%j.out
#SBATCH --error=debug_tmp/out/rebuild_serial_pymatching-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB

# build/pymatching predates the graph_cache format-version bump (GRAPH_CACHE_VERSION=2) --
# its serial baseline runs in benchmark_shmem_obs_cached.sh fail with "format version
# mismatch (file has 2, expected 1)" then abort (no --dem given as fallback). build/ is
# already configured (USE_SHMEM=OFF), so this just rebuilds the target -- no reconfigure.
source ~/.bash_profile
conda activate pymatching
cd ~/PyMatchingSHMEM
cmake --build build --target pymatching -j 16
