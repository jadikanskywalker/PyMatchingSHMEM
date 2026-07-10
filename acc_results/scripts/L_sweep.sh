#!/bin/bash
#SBATCH --job-name=acc_L_sweep
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=zen4
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=100GB

# Accuracy sweep for unit-checkpointed (bounded) extraction: --extraction_unit_size (L) and
# --extraction_buffer_levels (B). See plans/profiling-reveals-that-thread-buzzing-milner.md.
#
# Milestone 1 scope: single-observable DEM, ROUND division strategy, single PE. Multi-observable
# runs with the ROUND strategy are not fully tested in this codebase yet -- use a simple
# single-observable-patch DEM here, same as this script does, not a multi-obs surgery preset.
#
# L=0 is the baseline (checkpointing disabled, today's whole-root extraction); every other (L, B)
# is compared against it both for logical error rate and bit-exact equality. Because bounded
# extraction is an intentional, tunable accuracy/latency tradeoff (see Finding 1 in the plan doc),
# expect error_rate to converge exactly to the baseline once L/B are large enough relative to the
# code distance and error rate, with a small, smoothly graded degradation as L/B shrink -- not a
# step-function correctness break.

cd ~/PyMatchingSHMEM
source ~/.bash_profile
source ~/.bashrc
conda activate pymatching 2>/dev/null || true

set -eo pipefail

PROJECT_DIR=~/PyMatchingSHMEM
SHMEM_DECODER=$PROJECT_DIR/build_sos/pymatching
OSHRUN=${SWHOME}/sos_1.5_scalable/bin/oshrun
OUT=$PROJECT_DIR/acc_results/out

echo "Environment ready. OSHRUN=$OSHRUN"
ls "$SHMEM_DECODER" "$OSHRUN" || { echo "FATAL: binary not found"; exit 1; }

D=9
ROUNDS=200
P=0.005
SHOTS=1000
M=5
NTHREADS=16
TAG="1obs_d${D}_p${P//./}_r${ROUNDS}"

DEM_CACHE="$OUT/dems/${TAG}.dem"
DET="$OUT/events/${TAG}_${SHOTS}s.b8"
ACTUAL="$OUT/events/${TAG}_${SHOTS}s.01"

L_VALUES=(0 2 4 8 16 32)
B_VALUES=(0 1 2)

export OMP_NUM_THREADS=$NTHREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=true
export SHMEM_SYMMETRIC_SIZE=4G

mkdir -p "$OUT/dems" "$OUT/events" "$OUT/preds"

# Generate single-observable DEM + sample events (once, cached)
if [ ! -f "$DEM_CACHE" ] || [ ! -f "$DET" ]; then
    echo "=== Generating $TAG ==="
    "$OSHRUN" \
        -n 1 \
        --map-by ppr:1:node:PE=$NTHREADS \
        --bind-to core \
        "$SHMEM_DECODER" predict \
        --gen_code surface_code \
        --gen_task rotated_memory_x \
        --gen_distance $D \
        --gen_rounds $ROUNDS \
        --gen_num_obs 1 \
        --gen_depolarization $P \
        --dem_cache_path "$DEM_CACHE" \
        --gen_det_out "$DET" \
        --gen_obs_out "$ACTUAL" \
        --gen_sample_shots $SHOTS \
        --gen_sample_seed 42 \
        --in "$DET" \
        --in_format b8 \
        --out "$OUT/preds/pred_${TAG}_gencheck.01" \
        --out_format 01 \
        --rounds_per_partition $M \
        --task_division_strategy round \
        --use_threads \
        --extraction_unit_size 0
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

CSV=$OUT/L_sweep_results.csv
echo "L,B,matches,total,error_rate,bit_exact_vs_L0" > "$CSV"

BASELINE_PRED="$OUT/preds/pred_${TAG}_L0_B0_pe0.01"

for L in "${L_VALUES[@]}"; do
    for B in "${B_VALUES[@]}"; do
        # B is only meaningful when L>0; run it once (B=0) for the L=0 baseline row.
        if [ "$L" -eq 0 ] && [ "$B" -ne 0 ]; then
            continue
        fi

        # USE_SHMEM inserts "_pe<pid>" before the extension of --out (namespaced_main.cc), so the
        # actual file written is pred_..._peN.01, not pred_....01 -- account for that here.
        out_base="$OUT/preds/pred_${TAG}_L${L}_B${B}"
        pred="${out_base}_pe0.01"
        echo "=== L=$L B=$B ==="

        run_ok=true
        "$OSHRUN" \
            -n 1 \
            --map-by ppr:1:node:PE=$NTHREADS \
            --bind-to core \
            "$SHMEM_DECODER" predict \
            --dem "$DEM_CACHE" \
            --in "$DET" \
            --in_format b8 \
            --out "${out_base}.01" \
            --out_format 01 \
            --rounds_per_partition $M \
            --task_division_strategy round \
            --use_threads \
            --extraction_unit_size $L \
            --extraction_buffer_levels $B \
            || run_ok=false

        if $run_ok && [ -f "$pred" ]; then
            read matches total <<< $(compare_accuracy "$pred" "$ACTUAL")
            if [ "$total" -gt 0 ]; then
                error_rate=$(awk "BEGIN { printf \"%.6f\", 1 - $matches / $total }")
            else
                error_rate="NA"
            fi
            if [ -f "$BASELINE_PRED" ] && cmp -s "$pred" "$BASELINE_PRED"; then
                bit_exact="yes"
            else
                bit_exact="no"
            fi
        else
            echo "    FAILED (crash or error)"
            matches=0; total=0; error_rate="FAILED"; bit_exact="n/a"
        fi
        echo "  $matches/$total (error_rate=$error_rate, bit_exact_vs_L0=$bit_exact)"
        echo "$L,$B,$matches,$total,$error_rate,$bit_exact" >> "$CSV"
    done
done

echo "=== Done. Results in $CSV ==="
