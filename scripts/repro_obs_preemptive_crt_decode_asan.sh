#!/bin/bash
#SBATCH --job-name=repro_obs_preemptive_crt_decode_asan
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_out/repro_obs_preemptive_crt_decode_asan-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_out
export OMP_NUM_THREADS=4
export SHMEM_SYMMETRIC_SIZE=4G
export ASAN_OPTIONS=detect_leaks=0:halt_on_error=1

# ASan diagnostic for the crash found by verify_obs_preemptive_decode_regression_small.sh: preemptive
# OBS decode of real shots against the 4obs/2PE (all-CRT) DEM crashes with a SHMEM "not symmetric (0x1)"
# error inside wait_until_done, while the identical non-preemptive OBS run on the same shots succeeds.
# 1 shot only -- just need the crash to reproduce, not a full regression run.
DEM=~/PyMatchingSHMEM/testdems/error_model_4obs_d21_p001_231r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_4obs_d21_p001_231r_250s.b8

python3 - "$DET" "$DEM" verify_out/det_small1.b8 <<'EOF'
import sys
import stim
det_path, dem_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
dem = stim.DetectorErrorModel.from_file(dem_path)
num_dets = dem.num_detectors
data = stim.read_shot_data_file(path=det_path, format="b8", num_detectors=num_dets)
stim.write_shot_data_file(data=data[:1], path=out_path, format="b8", num_detectors=num_dets)
EOF

rm -rf out_parallel
oshrun -n 2 build_asan_throwaway/pymatching predict --dem "$DEM" \
  --in verify_out/det_small1.b8 --in_format b8 \
  --out verify_out/pred_obs_preemptive_decode_asan.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 1 \
  --extract_preemptively \
  > verify_out/repro_obs_preemptive_crt_decode_asan.stdout 2>&1
echo "exit=$?" >> verify_out/repro_obs_preemptive_crt_decode_asan.stdout
