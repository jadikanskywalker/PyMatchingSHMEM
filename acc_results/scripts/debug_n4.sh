#!/bin/bash
#SBATCH --job-name=debug_n4
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=zen4
#SBATCH --time=01:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
#SBATCH --mem=200GB

cd ~/PyMatchingSHMEM
source ~/.bash_profile
source ~/.bashrc
conda activate pymatching 2>/dev/null || true
export FI_VERBS_DEVICE_NAME="mlx5_2"

set -eo pipefail

PROJECT_DIR=~/PyMatchingSHMEM
DECODER=$PROJECT_DIR/build_sos_debug/pymatching
SHMEM_DECODER=$PROJECT_DIR/build_sos/pymatching
OSHRUN=${SWHOME}/sos_1.5_scalable/bin/oshrun
COMBINE_SCRIPT=$PROJECT_DIR/scripts/combine_results.py
OUT=$PROJECT_DIR/acc_results/out

echo "DECODER=$DECODER"
echo "OSHRUN=$OSHRUN"
ls "$DECODER" "$SHMEM_DECODER" "$OSHRUN" || { echo "FATAL: binary not found"; exit 1; }

D=21
ROUNDS=2058
P=0.001
SHOTS=100
M=21
NTHREADS=8
TAG="36obs_d21_p001_2058r"

DEM_CACHE="$OUT/dems/${TAG}.dem"
DET="$OUT/events/${TAG}_${SHOTS}s.b8"
ACTUAL="$OUT/events/${TAG}_${SHOTS}s.01"

export OMP_NUM_THREADS=$NTHREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export SHMEM_SYMMETRIC_SIZE=16G

# Generate 36obs DEM + sample events (once, cached)
if [ ! -f "$DEM_CACHE" ] || [ ! -f "$DET" ]; then
    echo "=== Generating $TAG via C++ gen ==="
    "$OSHRUN" \
        -n 1 \
        --map-by ppr:1:node:PE=$NTHREADS \
        --bind-to core \
        "$SHMEM_DECODER" predict \
        --gen_code surface_code \
        --gen_task rotated_memory_x \
        --gen_distance $D \
        --gen_rounds $ROUNDS \
        --gen_num_obs 36 \
        --gen_depolarization $P \
        --gen_surgery_preset 36obs \
        --dem_cache_path "$DEM_CACHE" \
        --gen_det_out "$DET" \
        --gen_obs_out "$ACTUAL" \
        --gen_sample_shots $SHOTS \
        --gen_sample_seed 42 \
        --in "$DET" \
        --in_format b8 \
        --out "$OUT/preds/pred_${TAG}_n1_shmem_k1.01" \
        --out_format 01 \
        --rounds_per_partition $M \
        --obs_coors_included \
        --cross_rank_fusion_window_size 1 \
        --task_division_strategy observable \
        --use_threads
    echo "  Cached: $DEM_CACHE"
else
    echo "=== Reusing cached $TAG ==="
fi

compare_accuracy() {
    local pred="$1" actual="$2"
    paste -d "" "$pred" "$actual" | awk '{
        n = length($0) / 2
        for (i = 1; i <= n; i++)
            if (substr($0, i, 1) == substr($0, n + i, 1)) c++
        t += n
    } END { printf "%d %d", c+0, t+0 }'
}

for k in 0 1; do
    echo "========================================"
    echo "=== n=4 k=$k 36obs (2 nodes, 2 ppn, debug) ==="
    echo "========================================"

    pred="$OUT/preds/debug_${TAG}_n4_k${k}.01"

    "$OSHRUN" \
        -n 4 \
        --map-by ppr:2:node:PE=$NTHREADS \
        --bind-to core \
        --report-bindings \
        "$DECODER" predict \
        --dem "$DEM_CACHE" \
        --in "$DET" \
        --in_format b8 \
        --out "$pred" \
        --out_format 01 \
        --rounds_per_partition $M \
        --obs_coors_included \
        --cross_rank_fusion_window_size $k \
        --task_division_strategy observable \
        --use_threads \
        || { echo "=== k=$k CRASHED ==="; continue; }

    python3 "$COMBINE_SCRIPT" "$pred" 4

    read matches total <<< $(compare_accuracy "$pred" "$ACTUAL")
    error_rate=$(awk "BEGIN { printf \"%.6f\", 1 - $matches / $total }")
    echo "=== k=$k result: $matches/$total (error_rate=$error_rate) ==="
done

echo "=== Done ==="
echo "Debug logs: out_parallel/p{0,1,2,3}/"
