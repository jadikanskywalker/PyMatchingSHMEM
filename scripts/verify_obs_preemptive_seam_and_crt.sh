#!/bin/bash
#SBATCH --job-name=verify_obs_preemptive_seam_and_crt
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_out/verify_obs_preemptive_seam_and_crt-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_out
export OMP_NUM_THREADS=4
export SHMEM_SYMMETRIC_SIZE=16G

# 9 observable patches / 3 PEs = 3 local patches per PE -- each PE gets local seams WITHIN its
# own 3 patches plus CrossRankTasks at the boundaries to neighboring PEs' patches. Exercises the
# "one local seam + one CRT on the same PE" case from now-its-time-to-jazzy-turing.md step 5.
DEM=~/PyMatchingSHMEM/testdems/error_model_9obs_d21_p001_672r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_9obs_d21_p001_672r_100s.b8

python3 - "$DET" "$DEM" verify_out/det_zero_9obs_3pe.b8 <<'EOF'
import sys
import stim
det_path, dem_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
dem = stim.DetectorErrorModel.from_file(dem_path)
num_dets = dem.num_detectors
data = stim.read_shot_data_file(path=det_path, format="b8", num_detectors=num_dets)
stim.write_shot_data_file(data=data[:0], path=out_path, format="b8", num_detectors=num_dets)
EOF

rm -rf out_parallel
oshrun -n 3 build_sos/pymatching predict --dem "$DEM" \
  --in verify_out/det_zero_9obs_3pe.b8 --in_format b8 \
  --out verify_out/pred_obs_preemptive_seam_and_crt.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 1 \
  --extract_preemptively \
  > verify_out/verify_obs_preemptive_seam_and_crt.stdout 2>&1
echo "exit=$?" >> verify_out/verify_obs_preemptive_seam_and_crt.stdout
