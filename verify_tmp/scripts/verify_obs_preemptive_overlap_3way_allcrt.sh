#!/bin/bash
#SBATCH --job-name=verify_obs_preemptive_overlap_3way_allcrt
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_tmp/out/verify_obs_preemptive_overlap_3way_allcrt-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_tmp/out
export OMP_NUM_THREADS=4
export SHMEM_SYMMETRIC_SIZE=16G

# Phase 1.6 test: 9obs star-topology DEM (all 16 seams touch observable 8), 9 PEs (1 observable each) --
# every seam touching observable 8 becomes a CrossRankTask for PE8. At L=12, seams 0 (oi=0,vb_right=3),
# 1 (oi=4,vb_right=7), and 2 (oi=1,vb_right=11) all land on unit_hi=0 (floor(3/12)=floor(7/12)=floor(11/12)=0)
# -- a genuine 3-way, all-CrossRankTask collision on PE8's own chain. Confirmed analytically during
# planning by bucketing seam_infos' vb_right by L; verify against the actual seam_infos DEBUG dump too.
# Expect success (exit=0); hand-verify via the DEBUG dump that PE8's chain shows three stacked
# CrossRankTasks (local_child -> CRT -> CRT -> CRT -> fusion) wrapped by exactly one Task.
DEM=~/PyMatchingSHMEM/testdems/error_model_9obs_d21_p001_672r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_9obs_d21_p001_672r_100s.b8

python3 - "$DET" "$DEM" verify_tmp/out/det_zero_3way_allcrt.b8 <<'EOF'
import sys
import stim
det_path, dem_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
dem = stim.DetectorErrorModel.from_file(dem_path)
num_dets = dem.num_detectors
data = stim.read_shot_data_file(path=det_path, format="b8", num_detectors=num_dets)
stim.write_shot_data_file(data=data[:0], path=out_path, format="b8", num_detectors=num_dets)
EOF

rm -rf out_parallel
oshrun -n 9 build_sos/pymatching predict --dem "$DEM" \
  --in verify_tmp/out/det_zero_3way_allcrt.b8 --in_format b8 \
  --out verify_tmp/out/pred_obs_preemptive_3way_allcrt.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 12 \
  --extract_preemptively \
  > verify_tmp/out/verify_obs_preemptive_3way_allcrt.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/verify_obs_preemptive_3way_allcrt.stdout
