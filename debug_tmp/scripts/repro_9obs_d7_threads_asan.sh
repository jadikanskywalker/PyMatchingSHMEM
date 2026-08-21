#!/bin/bash
#SBATCH --job-name=repro_9obs_d7_asan
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=64GB
#SBATCH --output=debug_tmp/out/repro_9obs_d7_asan-%j.out
#SBATCH --error=debug_tmp/out/repro_9obs_d7_asan-%j.err

source ~/.bash_profile
conda activate pymatching

export OMP_NUM_THREADS=9
export OMP_PLACES=cores
export OMP_PROC_BIND=true

DEM=~/PyMatchingSHMEM/testdems/error_model_9obs_d7_p01_337r.dem
DET=~/PyMatchingSHMEM/testdems/detection_events_9obs_d7_p01_337r_50s.b8

mkdir -p ~/PyMatchingSHMEM/debug_tmp/out/repro_9obs_d7_threads_asan
cd ~/PyMatchingSHMEM/debug_tmp/out/repro_9obs_d7_threads_asan

echo "Starting build_threads (ASan) repro for 9obs d7 337r, 9 threads, draw_frames..."
start=$(date +%s)
~/PyMatchingSHMEM/build_threads/pymatching predict \
    --dem "$DEM" \
    --in "$DET" \
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
