#!/bin/bash
#SBATCH --job-name=verify_obs_preemptive_overlap_2seams
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_out/verify_obs_preemptive_overlap_2seams-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_out
export OMP_NUM_THREADS=8

# Phase 1.5 test: at L=8, the 4obs DEM's seam1 (vb_right=9) and seam2 (vb_right=15) both land on
# unit_hi=1 (floor(9/8)=floor(15/8)=1) -- a genuine local-local same-unit_hi collision (both seams are
# between the same, single-PE-local obs0/obs1 pair). Previously this threw; now must succeed.
DEM=~/PyMatchingSHMEM/testdems/error_model_4obs_d21_p001_231r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_4obs_d21_p001_231r_250s.b8

python3 - "$DET" "$DEM" verify_out/det_zero_overlap2.b8 <<'EOF'
import sys
import stim
det_path, dem_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
dem = stim.DetectorErrorModel.from_file(dem_path)
num_dets = dem.num_detectors
data = stim.read_shot_data_file(path=det_path, format="b8", num_detectors=num_dets)
stim.write_shot_data_file(data=data[:0], path=out_path, format="b8", num_detectors=num_dets)
EOF

rm -rf out_parallel
oshrun -n 1 build_sos/pymatching predict --dem "$DEM" \
  --in verify_out/det_zero_overlap2.b8 --in_format b8 \
  --out verify_out/pred_obs_preemptive_overlap2.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 8 \
  --extract_preemptively \
  > verify_out/verify_obs_preemptive_overlap2.stdout 2>&1
echo "exit=$?" >> verify_out/verify_obs_preemptive_overlap2.stdout
