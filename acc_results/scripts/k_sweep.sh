#!/bin/bash
#SBATCH --job-name=acc_k_sweep
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
SERIAL_DECODER=~/PyMatching/build/pymatching
SHMEM_DECODER=$PROJECT_DIR/build_sos/pymatching
OSHRUN=${SWHOME}/sos_1.5_scalable/bin/oshrun
COMBINE_SCRIPT=$PROJECT_DIR/scripts/combine_results.py
OUT=$PROJECT_DIR/acc_results/out

echo "Environment ready. OSHRUN=$OSHRUN"
ls "$SHMEM_DECODER" "$SERIAL_DECODER" "$OSHRUN" || { echo "FATAL: binary not found"; exit 1; }

NUM_OBS=4
D=21
ROUNDS=231
P=0.001
SHOTS=250
M=21
NTHREADS=8
SPEC="0,1,0,21;0,1,63,21;0,1,126,21;0,1,189,21"
TAG="4obs_d21_p001_231r"

DEM_CACHE="$OUT/dems/${TAG}.dem"
DET="$OUT/events/${TAG}_${SHOTS}s.b8"
ACTUAL="$OUT/events/${TAG}_${SHOTS}s.01"

K_VALUES=(0 1 2)
N_PES_LIST=(1 2 4)

export OMP_NUM_THREADS=$NTHREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export SHMEM_SYMMETRIC_SIZE=4G

# Generate 4obs DEM + sample events (once, cached)
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
        --gen_num_obs $NUM_OBS \
        --gen_depolarization $P \
        --gen_surgery_spec "$SPEC" \
        --gen_surgery_duration $D \
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

CSV=$OUT/k_sweep_results.csv
echo "k,n_pes,matches,total,error_rate" > "$CSV"

for k in "${K_VALUES[@]}"; do
    for n in "${N_PES_LIST[@]}"; do
        pred="$OUT/preds/pred_${TAG}_n${n}_k${k}.01"

        echo "=== k=$k n=$n ==="

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
            if [ "$total" -gt 0 ]; then
                error_rate=$(awk "BEGIN { printf \"%.6f\", 1 - $matches / $total }")
            else
                error_rate="NA"
            fi
        else
            echo "    FAILED (segfault or error)"
            matches=0; total=0; error_rate="FAILED"
        fi
        echo "  $matches/$total (error_rate=$error_rate)"
        echo "$k,$n,$matches,$total,$error_rate" >> "$CSV"
    done
done

rm -f "$OUT/hostfile.txt"
echo "=== Done. Results in $CSV ==="
