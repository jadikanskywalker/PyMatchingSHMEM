#!/bin/bash
#SBATCH --job-name=verify_obs_extraction_2pe
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_tmp/out/verify_obs_extraction_2pe-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_tmp/out
export OMP_NUM_THREADS=2

DEM=~/PyMatchingSHMEM/testdems/error_model_4obs_d21_p001_231r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_4obs_d21_p001_231r_250s.b8

# Truncate to a handful of shots so the DEBUG=1 task-tree/trace output stays small and readable.
python3 - "$DET" "$DEM" verify_tmp/out/det_small_2pe.b8 <<'EOF'
import sys
import stim
det_path, dem_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
dem = stim.DetectorErrorModel.from_file(dem_path)
num_dets = dem.num_detectors
data = stim.read_shot_data_file(path=det_path, format="b8", num_detectors=num_dets)
n_keep = min(5, data.shape[0])
print(f"num_detectors={num_dets} total_shots={data.shape[0]} keeping={n_keep}")
stim.write_shot_data_file(data=data[:n_keep], path=out_path, format="b8", num_detectors=num_dets)
EOF

rm -rf out_parallel
export SHMEM_SYMMETRIC_SIZE=4G
oshrun -n 2 build_sos/pymatching predict --dem "$DEM" \
  --in verify_tmp/out/det_small_2pe.b8 --in_format b8 \
  --out verify_tmp/out/pred_obs_2pe.01 --out_format 01 \
  --task_division_strategy observable \
  --extraction_unit_size 2 \
  > verify_tmp/out/verify_obs_2pe.stdout 2>&1
echo "exit=$?" >> verify_tmp/out/verify_obs_2pe.stdout
