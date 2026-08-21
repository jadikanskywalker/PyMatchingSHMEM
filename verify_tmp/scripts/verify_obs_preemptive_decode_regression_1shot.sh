#!/bin/bash
#SBATCH --job-name=verify_obs_preemptive_decode_regression_small
#SBATCH --partition=zen4
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_tmp/out/verify_obs_preemptive_decode_regression_small-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_tmp/out
export OMP_NUM_THREADS=4
export SHMEM_SYMMETRIC_SIZE=4G

# Smaller/faster real-shot exercise of decode_shots()'s new Phase 2 preemptive-OBS logic than
# verify_obs_preemptive_decode_regression.sh (whose 9obs/d21/672r DEM proved too slow for a 20-min
# budget even at 100 shots -- non-preemptive ground truth alone didn't finish). Reuses the 4obs
# DEM/2PE topology from verify_obs_preemptive_overlap_crt_crt.sh (all seams become CRTs across the 2
# PEs), but with a small slice of REAL shots instead of 0, at extraction_unit_size=1 (finer chains,
# more deferred connectors exercised than that script's own extraction_unit_size=8).
DEM=~/PyMatchingSHMEM/testdems/error_model_4obs_d21_p001_231r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_4obs_d21_p001_231r_250s.b8

python3 - "$DET" "$DEM" verify_tmp/out/det_small1_v2.b8 <<'EOF'
import sys
import stim
det_path, dem_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
dem = stim.DetectorErrorModel.from_file(dem_path)
num_dets = dem.num_detectors
data = stim.read_shot_data_file(path=det_path, format="b8", num_detectors=num_dets)
stim.write_shot_data_file(data=data[:1], path=out_path, format="b8", num_detectors=num_dets)
EOF

rm -rf out_parallel

echo "Non-preemptive OBS (ground truth)..."
oshrun -n 2 build_sos/pymatching predict --dem "$DEM" \
  --in verify_tmp/out/det_small1_v2.b8 --in_format b8 \
  --out verify_tmp/out/pred_obs_nonpreemptive_1shot.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 1 \
  > verify_tmp/out/verify_obs_preemptive_decode_regression_1shot_nonpreemptive.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/verify_obs_preemptive_decode_regression_1shot_nonpreemptive.stdout

echo "Preemptive OBS (Phase 2 under test)..."
oshrun -n 2 build_sos/pymatching predict --dem "$DEM" \
  --in verify_tmp/out/det_small1_v2.b8 --in_format b8 \
  --out verify_tmp/out/pred_obs_preemptive_decode_1shot.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 1 \
  --extract_preemptively \
  > verify_tmp/out/verify_obs_preemptive_decode_regression_1shot_preemptive.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/verify_obs_preemptive_decode_regression_1shot_preemptive.stdout

echo "Diffing predictions..."
if diff -q verify_tmp/out/pred_obs_nonpreemptive_1shot.01 verify_tmp/out/pred_obs_preemptive_decode_1shot.01 > /dev/null; then
    echo "MATCH: preemptive OBS decode is bit-identical to non-preemptive OBS decode (20 shots)"
else
    echo "MISMATCH: preemptive OBS decode differs from non-preemptive OBS decode"
    diff verify_tmp/out/pred_obs_nonpreemptive_1shot.01 verify_tmp/out/pred_obs_preemptive_decode_1shot.01 | head -20
fi
