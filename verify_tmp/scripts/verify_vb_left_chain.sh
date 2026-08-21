#!/bin/bash
#SBATCH --job-name=verify_vb_left
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_tmp/out/verify_vb_left-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_tmp/out
export OMP_NUM_THREADS=8

rm -rf out_parallel
build_threads/pymatching predict --dem verify_tmp/out/tree.dem \
  --in verify_tmp/out/det.b8 --in_format b8 \
  --out verify_tmp/out/pred_vbleft_check.01 --out_format 01 \
  --extraction_unit_size 1 --extract_preemptively --draw_frames \
  > verify_tmp/out/vbleft_check.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/vbleft_check.stdout
