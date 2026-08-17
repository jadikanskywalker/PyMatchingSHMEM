#!/bin/bash
#SBATCH --job-name=repro_36obs_asan
#SBATCH --partition=zen4
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=512GB
#SBATCH --output=repro_36obs_asan-%j.out
#SBATCH --error=repro_36obs_asan-%j.err

# Args: [mode]  -- "cached" (default, uses --graph_cache_path, no --draw_frames since that's
# incompatible with a graph loaded from cache) or "uncached" (parses the DEM fresh every run,
# --draw_frames available since the DEM object itself is on hand).
mode=${1:-cached}

cd ~/PyMatchingSHMEM
source ~/.bash_profile
conda activate pymatching

export OMP_NUM_THREADS=64
export OMP_PLACES=cores
export OMP_PROC_BIND=true

DEM=testdems/error_model_36obs_d21_p001_2058r.dem
DET=testdems/detection_events_36obs_d21_p001_2058r_100s.b8
CACHE=testdems/graph_36obs_d21_p001_2058r.cache

outdir=repro_36obs_threads_asan_$mode
mkdir -p "$outdir"
cd "$outdir"

if [ "$mode" = "cached" ]; then
    extra_flags="--graph_cache_path ../$CACHE"
elif [ "$mode" = "uncached" ]; then
    extra_flags="--draw_frames"
else
    echo "Unknown mode: $mode (expected 'cached' or 'uncached')"
    exit 1
fi

echo "Starting build_threads (ASan) repro for 36obs d21 2058r, mode=$mode..."
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
    $extra_flags \
    > run.out 2> run.err
echo "exit code: $?"
end=$(date +%s)
echo "Completed in $((end-start)) seconds"
