#!/bin/bash
#SBATCH --job-name=gen_d_sweep
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=zen4
#SBATCH --time=05:00:00
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

# run_config <tag> <shots> <M> <K> <obs> <d> <p> <gen_flags...>
#   Step 1: oshrun -n 1 to generate DEM + sample events (cached for reuse)
#   Step 2: serial decode
#   Step 3: parallel decodes (n=2, 4)
run_config() {
    local tag="$1"; shift
    local shots="$1"; shift
    local M="$1"; shift
    local K="$1"; shift
    local obs="$1"; shift
    local d="$1"; shift
    local p="$1"; shift

    local cache_dem="$OUT/dems/${tag}.dem"
    local events="$OUT/events/${tag}_${shots}s.b8"
    local actual="$OUT/events/${tag}_${shots}s.01"

    export OMP_NUM_THREADS=$NTHREADS
    export OMP_PLACES=cores
    export OMP_PROC_BIND=true
    export SHMEM_SYMMETRIC_SIZE=16G

    # Step 1: Generate + sample (oshrun -n 1, decodes as side effect)
    if [ ! -f "$cache_dem" ] || [ ! -f "$events" ]; then
        echo "=== Generating $tag ==="
        "$OSHRUN" \
            -n 1 \
            --map-by ppr:1:node:PE=$NTHREADS \
            --bind-to core \
            "$SHMEM_DECODER" predict \
            "$@" \
            --dem_cache_path "$cache_dem" \
            --gen_det_out "$events" \
            --gen_obs_out "$actual" \
            --gen_sample_shots $shots \
            --gen_sample_seed 42 \
            --in "$events" \
            --in_format b8 \
            --out "$OUT/preds/pred_${tag}_n1_shmem_k${K}.01" \
            --out_format 01 \
            --rounds_per_partition $M \
            --obs_coors_included \
            --cross_rank_fusion_window_size $K \
            --task_division_strategy observable \
            --use_threads
        echo "  Cached DEM: $cache_dem"
    else
        echo "=== Reusing cached $tag ==="
    fi

    # Step 2: Serial baseline
    local pred="$OUT/preds/pred_${tag}_n1_serial.01"
    echo "  Decoding serial ..."
    "$SERIAL_DECODER" predict \
        --dem "$cache_dem" \
        --in "$events" \
        --in_format b8 \
        --out "$pred" \
        --out_format 01

    read matches total <<< $(compare_accuracy "$pred" "$actual")
    local error_rate=$(awk "BEGIN { printf \"%.6f\", 1 - $matches / $total }")
    echo "    serial: $matches/$total (error_rate=$error_rate)"
    echo "$obs,$d,$p,1,$matches,$total,$error_rate" >> "$CSV"

    # Step 3: Parallel decodes
    for n in 2 4; do
        pred="$OUT/preds/pred_${tag}_n${n}_k${K}.01"
        echo "  Decoding n=$n k=$K ..."

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
            --cross_rank_fusion_window_size $K \
            --task_division_strategy observable \
            --use_threads

        python3 "$COMBINE_SCRIPT" "$pred" $n

        read matches total <<< $(compare_accuracy "$pred" "$actual")
        error_rate=$(awk "BEGIN { printf \"%.6f\", 1 - $matches / $total }")
        echo "    n=$n: $matches/$total (error_rate=$error_rate)"
        echo "$obs,$d,$p,$n,$matches,$total,$error_rate" >> "$CSV"
    done
}

CSV=$OUT/gen_d_sweep_results.csv
echo "obs,d,p,n_pes,matches,total,error_rate" > "$CSV"

# ========================
# 4obs sweep (small d, 1000 shots) — C++ generator with surgery spec
# ========================
declare -A SPECS_4OBS
SPECS_4OBS[5]="0,1,0,5;0,1,15,5;0,1,30,5;0,1,45,5"
SPECS_4OBS[7]="0,1,0,7;0,1,21,7;0,1,42,7;0,1,63,7"
SPECS_4OBS[9]="0,1,0,9;0,1,27,9;0,1,54,9;0,1,81,9"

for d in 5 7 9; do
    rounds=$((11 * d))
    for p in 0.001 0.003 0.005 0.007 0.01; do
        p_label=$(echo "$p" | sed 's/0\.//')
        tag="4obs_d${d}_p${p_label}_${rounds}r"

        run_config "$tag" 1000 $d 1 4 $d $p \
            --gen_code surface_code \
            --gen_task rotated_memory_x \
            --gen_distance $d \
            --gen_rounds $rounds \
            --gen_num_obs 4 \
            --gen_depolarization $p \
            --gen_surgery_spec "${SPECS_4OBS[$d]}" \
            --gen_surgery_duration $d
    done
done

# ========================
# 36obs sweep (d=21, 100 shots) — C++ generator with preset
# ========================
for p in 0.001 0.003 0.005 0.007 0.01; do
    p_label=$(echo "$p" | sed 's/0\.//')
    tag="36obs_d21_p${p_label}_2058r"

    run_config "$tag" 100 21 1 36 21 $p \
        --gen_code surface_code \
        --gen_task rotated_memory_x \
        --gen_distance 21 \
        --gen_rounds 2058 \
        --gen_num_obs 36 \
        --gen_depolarization $p \
        --gen_surgery_preset 36obs
done

rm -f "$OUT/hostfile.txt"
echo "=== Done. Results in $CSV ==="
