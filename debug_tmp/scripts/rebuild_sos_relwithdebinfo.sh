#!/bin/bash
#SBATCH --job-name=rebuild_sos_relwithdebinfo
#SBATCH --output=debug_tmp/out/sbatch/rebuild_sos_relwithdebinfo-%j.out
#SBATCH --error=debug_tmp/out/sbatch/rebuild_sos_relwithdebinfo-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# rpc-96-1 runs a newer glibc (2.38) than the rest of the fleet (2.34) -- a build landing there
# links against __isoc23_strtol/fmodf@GLIBC_2.38 and then fails to run on any normal node.
#SBATCH --exclude=rpc-96-1

source ~/.bash_profile
conda activate pymatching

cd ~/PyMatchingSHMEM/build_sos
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo .
make clean
make pymatching -j 31
