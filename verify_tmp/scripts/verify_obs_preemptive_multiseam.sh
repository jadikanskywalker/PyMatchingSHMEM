#!/bin/bash
#SBATCH --job-name=verify_obs_preemptive_multiseam
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_tmp/out/verify_obs_preemptive_multiseam-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_tmp/out
export OMP_NUM_THREADS=8
export SHMEM_SYMMETRIC_SIZE=16G

# 9 observable patches, all local at 1 PE -- exercises multiple local seams touching DIFFERENT
# observable pairs simultaneously (unlike the 4obs DEM, whose 4 seams are all between the same
# obs0/obs1 pair). Construction-only (0 shots), same rationale as verify_obs_preemptive_construction.sh.
DEM=~/PyMatchingSHMEM/testdems/error_model_9obs_d21_p001_672r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_9obs_d21_p001_672r_100s.b8

python3 - "$DET" "$DEM" verify_tmp/out/det_zero_9obs.b8 <<'EOF'
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
  --in verify_tmp/out/det_zero_9obs.b8 --in_format b8 \
  --out verify_tmp/out/pred_obs_preemptive_multiseam.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 1 \
  --extract_preemptively \
  > verify_tmp/out/verify_obs_preemptive_multiseam.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/verify_obs_preemptive_multiseam.stdout
