#!/bin/bash
#SBATCH --job-name=verify_round_regression_phase16
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_out/verify_round_regression_phase16-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_out
export OMP_NUM_THREADS=8

# Phase 1.6 regression: real ROUND-partitioning preemptive decode (actual shots, not 0-shot
# construction-only), using build_sos (already rebuilt with the TaskBase*-widened Task type).
# Reuses the same circuit generation as verify_round_asan.sh but against build_sos, comparing
# predictions to ground truth to confirm ROUND's own preemptive chain (left_child/right_child now
# TaskBase*, but ROUND never attaches a CRT as a child) is completely unaffected.
if [ ! -f verify_out/round_circuit.dem ]; then
  stim gen \
    --rounds 15 \
    --distance 9 \
    --after_clifford_depolarization 0.05 \
    --code repetition_code \
    --task memory \
    > verify_out/round_circuit.stim
  stim analyze_errors \
    --decompose_errors \
    --fold_loops \
    --in verify_out/round_circuit.stim \
    > verify_out/round_circuit.dem
  stim detect \
    --in verify_out/round_circuit.stim \
    --shots 200 \
    --obs_out verify_out/round_actual_obs_flips.01 \
    --obs_out_format 01 \
    --out verify_out/round_det.b8 \
    --out_format b8
fi

rm -rf out_parallel
oshrun -n 1 build_sos/pymatching predict \
  --dem verify_out/round_circuit.dem \
  --in verify_out/round_det.b8 --in_format b8 \
  --out verify_out/round_pred_phase16.01 --out_format 01 \
  --extraction_unit_size 4 \
  --extract_preemptively \
  > verify_out/verify_round_regression_phase16.stdout 2>&1
echo "exit=$?" >> verify_out/verify_round_regression_phase16.stdout
diff -q verify_out/round_pred_phase16.01 verify_out/round_actual_obs_flips.01 >> verify_out/verify_round_regression_phase16.stdout 2>&1
echo "diff_rc=$?" >> verify_out/verify_round_regression_phase16.stdout
