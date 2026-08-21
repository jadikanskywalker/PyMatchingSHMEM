#!/bin/bash
#SBATCH --job-name=rebuild_sos
#SBATCH --output=debug_tmp/out/rebuild_sos-%j.out
#SBATCH --error=debug_tmp/out/rebuild_sos-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

source ~/.bash_profile
conda activate pymatching

cd ~/PyMatchingSHMEM/build_sos
make pymatching -j 31
