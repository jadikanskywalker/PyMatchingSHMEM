#!/bin/bash
#SBATCH --job-name=verify_obs_preemptive_overlap_local_crt
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_tmp/out/verify_obs_preemptive_overlap_local_crt-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_tmp/out
export OMP_NUM_THREADS=4
export SHMEM_SYMMETRIC_SIZE=16G

# Phase 1.5 test: 9obs DEM, 2 PEs (PE0 owns obs[0,5), PE1 owns obs[5,9)). At L=5, on PE1, seam index 3
# (oi=5, LOCAL to PE1 since 5 is on PE1) and seam index 4 (oi=2, a CrossRankTask for PE1 since 2 is on
# PE0) both end on unit_hi=3 (floor(15/5)=floor(19/5)=3) -- a genuine local-seam + CRT collision on the
# same PE, same observable (8)'s own chain. Previously blocked by the (now-relaxed) blanket assertion;
# now must succeed (no throw), and PE0 has no such collision at this L (sanity: confirms the assertion's
# narrowing is precise, not accidentally over- or under-permissive).
DEM=~/PyMatchingSHMEM/testdems/error_model_9obs_d21_p001_672r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_9obs_d21_p001_672r_100s.b8

python3 - "$DET" "$DEM" verify_tmp/out/det_zero_local_crt.b8 <<'EOF'
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
  --in verify_tmp/out/det_zero_local_crt.b8 --in_format b8 \
  --out verify_tmp/out/pred_obs_preemptive_local_crt.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 5 \
  --extract_preemptively \
  > verify_tmp/out/verify_obs_preemptive_local_crt.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/verify_obs_preemptive_local_crt.stdout
