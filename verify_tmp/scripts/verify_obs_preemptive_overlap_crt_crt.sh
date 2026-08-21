#!/bin/bash
#SBATCH --job-name=verify_obs_preemptive_overlap_crt_crt
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_tmp/out/verify_obs_preemptive_overlap_crt_crt-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_tmp/out
export OMP_NUM_THREADS=4
export SHMEM_SYMMETRIC_SIZE=4G

# Phase 1.6 test: 4obs DEM, 2 PEs (1 obs each) -- all 4 seams become CrossRankTasks for both PEs. At
# L=8, seam1 (vb_right=9) and seam2 (vb_right=15) both land on unit_hi=1 (floor(9/8)=floor(15/8)=1),
# and BOTH are CRTs -- this is exactly the case Phase 1.5 scoped out (CRT construction couldn't defer
# for a same-unit_hi sibling) and Phase 1.6 now supports (CRT construction stacks via cs.tip.node
# instead of committing immediately). Expect success (exit=0, no throw -- the narrowed assertion from
# Phase 1.5 was removed entirely in Phase 1.6).
DEM=~/PyMatchingSHMEM/testdems/error_model_4obs_d21_p001_231r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_4obs_d21_p001_231r_250s.b8

python3 - "$DET" "$DEM" verify_tmp/out/det_zero_crtcrt.b8 <<'EOF'
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
  --in verify_tmp/out/det_zero_crtcrt.b8 --in_format b8 \
  --out verify_tmp/out/pred_obs_preemptive_crtcrt.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 8 \
  --extract_preemptively \
  > verify_tmp/out/verify_obs_preemptive_crtcrt.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/verify_obs_preemptive_crtcrt.stdout
