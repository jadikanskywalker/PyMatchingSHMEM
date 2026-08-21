#!/bin/bash
#SBATCH --job-name=verify_obs_preemptive_overlap_3way_mixed
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_tmp/out/verify_obs_preemptive_overlap_3way_mixed-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_tmp/out
export OMP_NUM_THREADS=4
export SHMEM_SYMMETRIC_SIZE=16G

# Phase 1.6 test: 9obs DEM, 2 PEs (PE0 owns obs[0,5), PE1 owns obs[5,9), same split as
# verify_obs_preemptive_overlap_local_crt.sh). At L=10, on PE1, seam 2 (oi=1,vb_right=11,CRT), seam 3
# (oi=5,vb_right=15,LOCAL since 5 is on PE1), and seam 4 (oi=2,vb_right=19,CRT) all land on unit_hi=1
# (floor(11/10)=floor(15/10)=floor(19/10)=1) -- a genuine 3-way collision: 1 local seam + 2
# CrossRankTasks on the same PE, same observable (8)'s own chain. Confirmed analytically during
# planning by bucketing seam_infos' vb_right by L; verify against the actual seam_infos DEBUG dump too.
# Expect success (exit=0); hand-verify the local seam correctly has a CRT-stack operand on one side.
DEM=~/PyMatchingSHMEM/testdems/error_model_9obs_d21_p001_672r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_9obs_d21_p001_672r_100s.b8

python3 - "$DET" "$DEM" verify_tmp/out/det_zero_3way_mixed.b8 <<'EOF'
import sys
import stim
det_path, dem_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
dem = stim.DetectorErrorModel.from_file(dem_path)
num_dets = dem.num_detectors
data = stim.read_shot_data_file(path=det_path, format="b8", num_detectors=num_dets)
stim.write_shot_data_file(data=data[:0], path=out_path, format="b8", num_detectors=num_dets)
EOF

rm -rf out_parallel
oshrun -n 2 build_sos/pymatching predict --dem "$DEM" \
  --in verify_tmp/out/det_zero_3way_mixed.b8 --in_format b8 \
  --out verify_tmp/out/pred_obs_preemptive_3way_mixed.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 10 \
  --extract_preemptively \
  > verify_tmp/out/verify_obs_preemptive_3way_mixed.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/verify_obs_preemptive_3way_mixed.stdout
