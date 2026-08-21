#!/bin/bash
#SBATCH --job-name=check_partitions_3pe_asan
#SBATCH --partition=zen4
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_tmp/out/check_partitions_3pe_asan-%j.out
source ~/.bash_profile
conda activate pymatching
cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
export OMP_NUM_THREADS=2
export ASAN_OPTIONS=detect_leaks=0:halt_on_error=1
rm -rf out_parallel
oshrun -n 3 build_asan_throwaway/pymatching predict \
  --dem verify_tmp/out/round_circuit.dem \
  --in verify_tmp/out/round_det.b8 --in_format b8 \
  --out verify_tmp/out/round_pred_3pe_asan.01 --out_format 01 \
  --rounds_per_partition 2 \
  --extraction_unit_size 4 \
  --extract_preemptively \
  > verify_tmp/out/check_partitions_3pe_asan.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/check_partitions_3pe_asan.stdout
