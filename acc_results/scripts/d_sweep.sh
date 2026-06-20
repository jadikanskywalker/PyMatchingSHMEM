#!/bin/bash
#SBATCH --job-name=acc_d_sweep
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=zen4
#SBATCH --time=03:00:00
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

NTHREADS=8
K=1

export OMP_NUM_THREADS=$NTHREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export SHMEM_SYMMETRIC_SIZE=4G

declare -A SURGERY_SPECS
SURGERY_SPECS[5]="0,1,0,5;0,1,15,5;0,1,30,5;0,1,45,5"
SURGERY_SPECS[7]="0,1,0,7;0,1,21,7;0,1,42,7;0,1,63,7"
SURGERY_SPECS[9]="0,1,0,9;0,1,27,9;0,1,54,9;0,1,81,9"

compare_accuracy() {
    local pred="$1" actual="$2"
    paste -d "" "$pred" "$actual" | awk '{
        n = length($0) / 2
        for (i = 1; i <= n; i++)
            if (substr($0, i, 1) == substr($0, n + i, 1)) c++
        t += n
    } END { printf "%d %d", c+0, t+0 }'
}

CSV=$OUT/d_sweep_results.csv
echo "d,p,n_pes,matches,total,error_rate" > "$CSV"

for d in 5 7 9; do
    rounds=$((11 * d))
    M=$d
    spec="${SURGERY_SPECS[$d]}"
    shots=1000

    for p in 0.001 0.003 0.005 0.007 0.01; do
        p_label=$(echo "$p" | sed 's/0\.//')
        tag="4obs_d${d}_p${p_label}_${rounds}r"

        cache_dem="$OUT/dems/${tag}.dem"
        det="$OUT/events/${tag}_${shots}s.b8"
        actual="$OUT/events/${tag}_${shots}s.01"

        # Generate + sample once via C++ gen
        if [ ! -f "$cache_dem" ] || [ ! -f "$det" ]; then
            echo "=== Generating $tag ==="
            "$OSHRUN" \
                -n 1 \
                --map-by ppr:1:node:PE=$NTHREADS \
                --bind-to core \
                "$SHMEM_DECODER" predict \
                --gen_code surface_code \
                --gen_task rotated_memory_x \
                --gen_distance $d \
                --gen_rounds $rounds \
                --gen_num_obs 4 \
                --gen_depolarization $p \
                --gen_surgery_spec "$spec" \
                --gen_surgery_duration $d \
                --dem_cache_path "$cache_dem" \
                --gen_det_out "$det" \
                --gen_obs_out "$actual" \
                --gen_sample_shots $shots \
                --gen_sample_seed 42 \
                --in "$det" \
                --in_format b8 \
                --out "$OUT/preds/pred_${tag}_n1_shmem_k${K}.01" \
                --out_format 01 \
                --rounds_per_partition $M \
                --obs_coors_included \
                --cross_rank_fusion_window_size $K \
                --task_division_strategy observable \
                --use_threads
        else
            echo "=== Reusing cached $tag ==="
        fi

        # Decode: serial, n=2, n=4
        for n in 1 2 4; do
            pred="$OUT/preds/pred_${tag}_n${n}_k${K}.01"
            echo "  Decoding n=$n k=$K ..."

            run_ok=true
            if [ "$n" -eq 1 ]; then
                "$SERIAL_DECODER" predict \
                    --dem "$cache_dem" \
                    --in "$det" \
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
                    --dem "$cache_dem" \
                    --in "$det" \
                    --in_format b8 \
                    --out "$pred" \
                    --out_format 01 \
                    --rounds_per_partition $M \
                    --obs_coors_included \
                    --cross_rank_fusion_window_size $K \
                    --task_division_strategy observable \
                    --use_threads \
                    || run_ok=false

                if $run_ok; then
                    python3 "$COMBINE_SCRIPT" "$pred" $n
                fi
            fi

            if $run_ok && [ -f "$pred" ]; then
                read matches total <<< $(compare_accuracy "$pred" "$actual")
                if [ "$total" -gt 0 ]; then
                    error_rate=$(awk "BEGIN { printf \"%.6f\", 1 - $matches / $total }")
                else
                    error_rate="NA"
                fi
            else
                echo "    FAILED"
                matches=0; total=0; error_rate="FAILED"
            fi
            echo "    $matches/$total (error_rate=$error_rate)"
            echo "$d,$p,$n,$matches,$total,$error_rate" >> "$CSV"
        done
    done
done

echo "=== Done. Results in $CSV ==="
