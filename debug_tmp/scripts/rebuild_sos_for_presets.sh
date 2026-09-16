#!/bin/bash
#SBATCH --job-name=rebuild_sos_for_presets
#SBATCH --output=debug_tmp/out/sbatch/rebuild_sos_for_presets-%j.out
#SBATCH --error=debug_tmp/out/sbatch/rebuild_sos_for_presets-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB

# One-off rebuild before the gen_{72,128,256,144}obs.sh jobs, so they don't race on the same
# build_sos directory if run concurrently. Picks up the new preset_144obs/128obs/256obs code.
source ~/.bash_profile
conda activate pymatching
cd ~/PyMatchingSHMEM
cmake --build build_sos --target pymatching -j 16
