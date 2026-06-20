#!/bin/bash
#SBATCH --job-name=gen_k_sweep
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=zen4
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
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

NTHREADS=8
SHOTS=100
M=21

# Hostfile
if command -v srun &>/dev/null && [ -n "${SLURM_JOB_ID:-}" ]; then
    hosts=$(srun hostname | sort | uniq | paste -sd, -)
else
    hosts=$(hostname)
fi
hostfile=$OUT/hostfile.txt
> "$hostfile"
for host in ${hosts//,/ }; do
    echo "$host slots=4" >> "$hostfile"
done

compare_accuracy() {
    local pred="$1" actual="$2"
    paste -d "" "$pred" "$actual" | awk '{
        n = length($0) / 2
        for (i = 1; i <= n; i++)
            if (substr($0, i, 1) == substr($0, n + i, 1)) c++
        t += n
    } END { printf "%d %d", c+0, t+0 }'
}

CSV=$OUT/gen_k_sweep_results.csv
echo "obs,d,p,k,n_pes,matches,total,error_rate" > "$CSV"

K_VALUES=(0 1 2)
N_PES_LIST=(1 2 4)

# For each p, generate DEM + events once, then sweep k × n_pes
for p in 0.01; do
    p_label=$(echo "$p" | sed 's/0\.//')
    tag="36obs_d21_p${p_label}_2058r"

    cache_dem="$OUT/dems/${tag}.dem"
    events="$OUT/events/${tag}_${SHOTS}s.b8"
    actual="$OUT/events/${tag}_${SHOTS}s.01"

    export OMP_NUM_THREADS=$NTHREADS
    export OMP_PLACES=cores
    export OMP_PROC_BIND=true
    export SHMEM_SYMMETRIC_SIZE=16G

    # Generate + sample once (reuse if cached)
    if [ ! -f "$cache_dem" ] || [ ! -f "$events" ]; then
        echo "=== Generating $tag ==="
        "$OSHRUN" \
            -n 1 \
            --map-by ppr:1:node:PE=$NTHREADS \
            --bind-to core \
            "$SHMEM_DECODER" predict \
            --gen_code surface_code \
            --gen_task rotated_memory_x \
            --gen_distance 21 \
            --gen_rounds 2058 \
            --gen_num_obs 36 \
            --gen_depolarization $p \
            --gen_surgery_preset 36obs \
            --dem_cache_path "$cache_dem" \
            --gen_det_out "$events" \
            --gen_obs_out "$actual" \
            --gen_sample_shots $SHOTS \
            --gen_sample_seed 42 \
            --in "$events" \
            --in_format b8 \
            --out "$OUT/preds/pred_${tag}_n1_shmem_k1.01" \
            --out_format 01 \
            --rounds_per_partition $M \
            --obs_coors_included \
            --cross_rank_fusion_window_size 1 \
            --task_division_strategy observable \
            --use_threads
        echo "  Cached: $cache_dem"
    else
        echo "=== Reusing cached $tag ==="
    fi

    # Sweep k × n_pes
    for k in "${K_VALUES[@]}"; do
        for n in "${N_PES_LIST[@]}"; do
            pred="$OUT/preds/pred_${tag}_n${n}_k${k}.01"

            echo "  p=$p k=$k n=$n ..."

            if [ "$n" -eq 1 ]; then
                "$SERIAL_DECODER" predict \
                    --dem "$cache_dem" \
                    --in "$events" \
                    --in_format b8 \
                    --out "$pred" \
                    --out_format 01
            else
                "$OSHRUN" \
                    -n $n \
                    --map-by ppr:${n}:node:PE=$NTHREADS \
                    --bind-to core \
                    "$SHMEM_DECODER" predict \
                    --dem "$cache_dem" \
                    --in "$events" \
                    --in_format b8 \
                    --out "$pred" \
                    --out_format 01 \
                    --rounds_per_partition $M \
                    --obs_coors_included \
                    --cross_rank_fusion_window_size $k \
                    --task_division_strategy observable \
                    --use_threads

                python3 "$COMBINE_SCRIPT" "$pred" $n
            fi

            read matches total <<< $(compare_accuracy "$pred" "$actual")
            if [ "$total" -gt 0 ]; then
                error_rate=$(awk "BEGIN { printf \"%.6f\", 1 - $matches / $total }")
            else
                error_rate="NA"
            fi
            echo "    $matches/$total (error_rate=$error_rate)"
            echo "36,21,$p,$k,$n,$matches,$total,$error_rate" >> "$CSV"
        done
    done
done

rm -f "$OUT/hostfile.txt"
echo "=== Done. Results in $CSV ==="
