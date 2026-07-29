#!/bin/bash
#SBATCH --job-name=verify_obs_preemptive_overlap_cascade
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_out/verify_obs_preemptive_overlap_cascade-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_out
export OMP_NUM_THREADS=8

# Phase 1.5 stress test: the 9obs DEM's topology is a star (all 16 seams have oj=8). At L=8, seven
# DIFFERENT pairs of these seams collide on the same unit_hi simultaneously (all local-local, since
# 1 PE owns all 9 patches) -- observable 8's own chain must correctly cascade through all of them,
# routing each seam's right_obs_parent independently. Previously this threw; now must succeed.
DEM=~/PyMatchingSHMEM/testdems/error_model_9obs_d21_p001_672r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_9obs_d21_p001_672r_100s.b8

python3 - "$DET" "$DEM" verify_out/det_zero_cascade.b8 <<'EOF'
import sys
import stim
det_path, dem_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
dem = stim.DetectorErrorModel.from_file(dem_path)
num_dets = dem.num_detectors
data = stim.read_shot_data_file(path=det_path, format="b8", num_detectors=num_dets)
stim.write_shot_data_file(data=data[:0], path=out_path, format="b8", num_detectors=num_dets)
EOF

rm -rf out_parallel
export SHMEM_SYMMETRIC_SIZE=16G
oshrun -n 1 build_sos/pymatching predict --dem "$DEM" \
  --in verify_out/det_zero_cascade.b8 --in_format b8 \
  --out verify_out/pred_obs_preemptive_cascade.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 8 \
  --extract_preemptively \
  > verify_out/verify_obs_preemptive_cascade.stdout 2>&1
echo "exit=$?" >> verify_out/verify_obs_preemptive_cascade.stdout
