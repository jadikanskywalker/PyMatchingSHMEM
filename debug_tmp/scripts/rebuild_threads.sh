#!/bin/bash
#SBATCH --job-name=rebuild_threads
#SBATCH --output=debug_tmp/out/sbatch/rebuild_threads-%j.out
#SBATCH --error=debug_tmp/out/sbatch/rebuild_threads-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB

source ~/.bash_profile
conda activate pymatching
cd ~/PyMatchingSHMEM
cmake --build build_threads --target pymatching -j 16
