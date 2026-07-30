#!/bin/bash
#SBATCH --job-name=repro_4obs_nullptr_parents
#SBATCH --partition=zen4
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/verify_out/repro_4obs_nullptr_parents-%j.out

source ~/.bash_profile
conda activate pymatching

cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM
mkdir -p verify_out
export OMP_NUM_THREADS=2
export SHMEM_SYMMETRIC_SIZE=16G

# Repro for user's testrun_shmem.sh report: 4obs, 63 rounds, M=16 (rounds_per_partition), 4 PEs,
# L=4, --extract_preemptively --task_division_strategy observable -- 3 seams reported (0-2, 0-1, 0-2)
# but root fusions for obs0/1/2 show parent==nullptr and no CRTs appear. Same gen_multi_obs.py
# generator as testrun_shmem.sh, minus --draw_frames/--use_threads to isolate the construction path.
python3 scripts/gen_multi_obs.py \
    --num_observables 4 \
    --rounds 63 \
    --distance 7 \
    --after_clifford_depolarization 0.1 \
    --code repetition_code \
    --task memory \
    --num_surgery_gates 3 \
    --surgery_duration 5 \
    --circuit_out verify_out/repro4obs_circuit.stim \
    > verify_out/repro4obs_error_model.dem

python3 - verify_out/repro4obs_error_model.dem verify_out/repro4obs_det_zero.b8 <<'EOF'
import sys
import numpy as np
import stim
dem_path, out_path = sys.argv[1], sys.argv[2]
dem = stim.DetectorErrorModel.from_file(dem_path)
num_dets = dem.num_detectors
stim.write_shot_data_file(data=np.zeros((0, num_dets), dtype=bool), path=out_path,
                           format="b8", num_detectors=num_dets)
EOF

rm -rf out_parallel
oshrun -n 4 build_sos/pymatching predict --dem verify_out/repro4obs_error_model.dem \
  --in verify_out/repro4obs_det_zero.b8 --in_format b8 \
  --out verify_out/repro4obs_pred.01 --out_format 01 \
  --rounds_per_partition 16 \
  --obs_coors_included \
  --cross_rank_fusion_window_size 1 \
  --task_division_strategy observable \
  --extraction_unit_size 4 \
  --extract_preemptively \
  > verify_out/repro_4obs_nullptr_parents.stdout 2>&1
echo "exit=$?" >> verify_out/repro_4obs_nullptr_parents.stdout
