#!/bin/bash
#SBATCH --job-name=verify_obs_preemptive_decode_regression
#SBATCH --partition=zen4
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_tmp/out/verify_obs_preemptive_decode_regression-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_tmp/out
export OMP_NUM_THREADS=4
export SHMEM_SYMMETRIC_SIZE=16G

# First real-shot (non-zero) exercise of decode_shots()'s new Phase 2 preemptive-OBS logic
# (now-its-time-to-jazzy-turing.md): fixed-slot seam routing, Bug A finalize-skip, Bug B CRT-chain
# decoupling. Construction-only DEBUG=0-shot dumps (verify_obs_preemptive_*.sh) never exercised any
# of this -- everything here only matters once real shots run through process_timeline_until_completion.
#
# 9obs/3PE/extraction_unit_size=1 mirrors verify_obs_preemptive_seam_and_crt.sh's topology (each PE
# gets local seams within its own 3 patches plus CRTs to neighboring PEs) -- exercises both local-seam
# and CRT decode paths on the same run. Ground truth: the same 100 shots decoded via non-preemptive OBS.
DEM=~/PyMatchingSHMEM/testdems/error_model_9obs_d21_p001_672r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_9obs_d21_p001_672r_100s.b8

rm -rf out_parallel

echo "Non-preemptive OBS (ground truth)..."
oshrun -n 3 build_sos/pymatching predict --dem "$DEM" \
  --in "$DET" --in_format b8 \
  --out verify_tmp/out/pred_obs_nonpreemptive.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 1 \
  > verify_tmp/out/verify_obs_preemptive_decode_regression_nonpreemptive.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/verify_obs_preemptive_decode_regression_nonpreemptive.stdout

echo "Preemptive OBS (Phase 2 under test)..."
oshrun -n 3 build_sos/pymatching predict --dem "$DEM" \
  --in "$DET" --in_format b8 \
  --out verify_tmp/out/pred_obs_preemptive_decode.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 1 \
  --extract_preemptively \
  > verify_tmp/out/verify_obs_preemptive_decode_regression_preemptive.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/verify_obs_preemptive_decode_regression_preemptive.stdout

echo "Diffing predictions..."
if diff -q verify_tmp/out/pred_obs_nonpreemptive.01 verify_tmp/out/pred_obs_preemptive_decode.01 > /dev/null; then
    echo "MATCH: preemptive OBS decode is bit-identical to non-preemptive OBS decode (100 shots)"
else
    echo "MISMATCH: preemptive OBS decode differs from non-preemptive OBS decode"
    diff verify_tmp/out/pred_obs_nonpreemptive.01 verify_tmp/out/pred_obs_preemptive_decode.01 | head -20
fi
