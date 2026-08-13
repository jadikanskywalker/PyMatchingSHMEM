#!/bin/bash
#SBATCH --job-name=repro_36obs_asan
#SBATCH --partition=zen4
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=256GB
#SBATCH --output=repro_36obs_asan-%j.out
#SBATCH --error=repro_36obs_asan-%j.err

cd ~/PyMatchingSHMEM
source ~/.bash_profile
conda activate pymatching

export OMP_NUM_THREADS=64
export OMP_PLACES=cores
export OMP_PROC_BIND=true

DEM=testdems/error_model_36obs_d21_p001_2058r.dem
DET=testdems/detection_events_36obs_d21_p001_2058r_100s.b8

mkdir -p repro_36obs_threads
cd repro_36obs_threads

echo "Starting build_threads (ASan) repro..."
start=$(date +%s)
../build_threads/pymatching predict \
    --dem ../$DEM \
    --in ../$DET \
    --in_format b8 \
    --out predicted.01 \
    --out_format 01 \
    --rounds_per_partition 21 \
    --seam_buffer_size 1 \
    --task_division_strategy observable \
    --extraction_unit_size 4 \
    --extract_preemptively \
    --use_threads \
    --draw_frames \
    > run.out 2> run.err
echo "exit code: $?"
end=$(date +%s)
echo "Completed in $((end-start)) seconds"
