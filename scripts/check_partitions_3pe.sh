#!/bin/bash
#SBATCH --job-name=check_partitions_3pe
#SBATCH --partition=zen4
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/mnt/DISCL/home/jadhicks/PyMatchingSHMEM/run/check_partitions_3pe-%j.out

source ~/.bashrc
conda activate pymatching
cd /mnt/DISCL/home/jadhicks/PyMatchingSHMEM/run
export OMP_NUM_THREADS=1
rm -rf out_parallel

stim gen \
    --rounds 128 \
    --distance 5 \
    --after_clifford_depolarization 0.1 \
    --code repetition_code \
    --task memory \
    > circuit.stim
stim analyze_errors \
    --decompose_errors \
    --fold_loops \
    --in circuit.stim \
    > error_model.dem
stim detect \
    --in circuit.stim \
    --shots 10 \
    --obs_out actual_obs_flips.01 \
    --obs_out_format 01 \
    --out detection_events.b8 \
    --out_format b8

oshrun -n 3 ../build_sos/pymatching predict \
  --dem error_model.dem \
  --in detection_events.b8 --in_format b8 \
  --out predicted_obs_flips.01 --out_format 01 \
  --rounds_per_partition 8 \
  --extraction_unit_size 2 \
  --extract_preemptively \
  > check_partitions_3pe.stdout 2>&1
echo "exit=$?" >> check_partitions_3pe.stdout
