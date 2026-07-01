#!/bin/bash
#SBATCH --job-name=64obs_acc
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=zen4
#SBATCH --time=36:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
#SBATCH --mem=750GB

# Usage: sbatch run_64obs.sh <error_probability>
# Example: sbatch run_64obs.sh 0.005

cd ~/PyMatchingSHMEM
source ~/.bash_profile
source ~/.bashrc
conda activate pymatching 2>/dev/null || true
export FI_VERBS_DEVICE_NAME="mlx5_2"

set -eo pipefail

P="${1:?Usage: sbatch $0 <error_probability>}"
p_label=$(echo "$P" | sed 's/0\.//')

PROJECT_DIR=~/PyMatchingSHMEM
SERIAL_DECODER=~/PyMatching/build/pymatching
SHMEM_DECODER=$PROJECT_DIR/build_sos/pymatching
OSHRUN=${SWHOME}/sos_1.5_scalable/bin/oshrun
COMBINE_SCRIPT=$PROJECT_DIR/scripts/combine_results.py
OUT=$PROJECT_DIR/acc_results/out

SHOTS=1000
M=21
NTHREADS=8
TAG="64obs_d21_p${p_label}_704r"

DEM_CACHE="$OUT/dems/${TAG}.dem"
DET="$OUT/events/${TAG}_${SHOTS}s.b8"
ACTUAL="$OUT/events/${TAG}_${SHOTS}s.01"

echo "=== 64obs p=$P ($TAG) ==="
echo "OSHRUN=$OSHRUN"
ls "$SHMEM_DECODER" "$SERIAL_DECODER" "$OSHRUN" || { echo "FATAL: binary not found"; exit 1; }

export OMP_NUM_THREADS=$NTHREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export SHMEM_SYMMETRIC_SIZE=16G

# Generate DEM (if not cached) + sample 1000 shots
if [ ! -f "$DEM_CACHE" ] || [ ! -f "$DET" ]; then
    echo "=== Generating + sampling $TAG ==="
    "$OSHRUN" \
        -n 1 \
        --map-by ppr:1:node:PE=$NTHREADS \
        --bind-to core \
        "$SHMEM_DECODER" predict \
        --gen_code surface_code \
        --gen_task rotated_memory_x \
        --gen_distance 21 \
        --gen_rounds 704 \
        --gen_num_obs 64 \
        --gen_depolarization $P \
        --gen_surgery_preset 64obs \
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
    echo "=== Reusing cached DEM + events ==="
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

CSV=$OUT/run_64obs_p${p_label}_results.csv
echo "obs,d,p,k,n_pes,matches,total,error_rate" > "$CSV"

for k in 0 1 2; do
    for n in 1 2 4; do
        pred="$OUT/preds/pred_${TAG}_n${n}_k${k}.01"
        echo "=== p=$P k=$k n=$n ==="

        run_ok=true
        if [ "$n" -eq 1 ]; then
            "$SERIAL_DECODER" predict \
                --dem "$DEM_CACHE" \
                --in "$DET" \
                --in_format b8 \
                --out "$pred" \
                --out_format 01 \
                || run_ok=false
        else
            "$OSHRUN" \
                -n $n \
                --map-by ppr:2:node:PE=$NTHREADS \
                --bind-to core \
                "$SHMEM_DECODER" predict \
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
                || run_ok=false

            if $run_ok; then
                python3 "$COMBINE_SCRIPT" "$pred" $n
            fi
        fi

        if $run_ok && [ -f "$pred" ]; then
            read matches total <<< $(compare_accuracy "$pred" "$ACTUAL")
            error_rate=$(awk "BEGIN { printf \"%.6f\", 1 - $matches / $total }")
        else
            echo "    FAILED"
            matches=0; total=0; error_rate="FAILED"
        fi
        echo "  $matches/$total (error_rate=$error_rate)"
        echo "64,21,$P,$k,$n,$matches,$total,$error_rate" >> "$CSV"
    done
done

echo "=== Done. Results in $CSV ==="
