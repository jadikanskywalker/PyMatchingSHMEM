#!/bin/bash
#SBATCH --job-name=verify_round_asan
#SBATCH --partition=zen4
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_out/verify_round_asan-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_out
export OMP_NUM_THREADS=8
export ASAN_OPTIONS=detect_leaks=0:halt_on_error=1

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

FAILS=0
for i in $(seq 1 20); do
  rm -rf out_parallel
  oshrun -n 1 build_asan_throwaway/pymatching predict \
    --dem verify_out/round_circuit.dem \
    --in verify_out/round_det.b8 --in_format b8 \
    --out verify_out/round_pred_${i}.01 --out_format 01 \
    --extraction_unit_size 4 \
    --extract_preemptively \
    > verify_out/round_asan_run_${i}.stdout 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "RUN $i FAILED rc=$rc"
    FAILS=$((FAILS+1))
  else
    echo "RUN $i ok"
    diff -q verify_out/round_pred_${i}.01 verify_out/round_actual_obs_flips.01 >> verify_out/round_asan_diff.log 2>&1
  fi
done
echo "TOTAL_FAILS=$FAILS" | tee -a verify_out/round_asan_summary.txt
