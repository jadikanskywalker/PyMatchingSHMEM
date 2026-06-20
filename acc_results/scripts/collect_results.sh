#!/bin/bash
# Collect and display accuracy results.
# Usage:
#   ./collect_results.sh                          # default: searches testdems/ and acc_results/out/events/
#   ./collect_results.sh --gt-dir path/to/flips   # specify ground truth folder

set -euo pipefail

PROJECT_DIR=~/PyMatchingSHMEM
OUT=$PROJECT_DIR/acc_results/out

# Parse options
GT_DIRS=("$PROJECT_DIR/testdems" "$OUT/events")
while [[ $# -gt 0 ]]; do
    case "$1" in
        --gt-dir)
            GT_DIRS=("$2")
            shift 2
            ;;
        *)
            echo "Usage: $0 [--gt-dir <path>]" >&2
            exit 1
            ;;
    esac
done

# --- CSV results ---
echo "===== Sweep Result CSVs ====="
for csv in "$OUT"/*_results.csv; do
    if [ -f "$csv" ]; then
        echo
        echo "--- $(basename "$csv") ---"
        column -t -s',' "$csv"
    fi
done

# --- Per-file accuracy against ground truth ---
echo
echo "===== Per-File Accuracy (preds vs ground truth) ====="

compare_accuracy() {
    local pred="$1" gt="$2"
    paste -d "" "$pred" "$gt" | awk '{
        n = length($0) / 2
        for (i = 1; i <= n; i++)
            if (substr($0, i, 1) == substr($0, n + i, 1)) c++
        t += n
    } END { printf "%d %d", c+0, t+0 }'
}

# Find ground truth for a prediction file by matching the tag in its name.
# Searches GT_DIRS for actual_obs_flips_<tag>_<shots>s.01 or <tag>_<shots>s.01
find_gt() {
    local pred_name="$1"
    # Extract tag: pred_<tag>_n<N>_k<K>.01 -> <tag>
    local tag=$(echo "$pred_name" | sed 's/^pred_//; s/_n[0-9]*_k[0-9]*\.01$//')

    for dir in "${GT_DIRS[@]}"; do
        # acc_results/out/events style: <tag>_<shots>s.01
        for gt in "$dir"/${tag}_*s.01; do
            if [ -f "$gt" ]; then
                echo "$gt"
                return
            fi
        done
        # testdems style: actual_obs_flips_<tag>.01
        for gt in "$dir"/actual_obs_flips_${tag}*.01; do
            if [ -f "$gt" ]; then
                echo "$gt"
                return
            fi
        done
    done
}

printf "%-60s %s\n" "FILE" "MATCH/TOTAL"
printf "%0.s-" {1..75}; echo

for pred in "$OUT"/preds/pred_*.01; do
    [ -f "$pred" ] || continue
    name=$(basename "$pred")

    # Skip per-PE files
    [[ "$name" == *_pe[0-9]* ]] && continue

    gt=$(find_gt "$name")
    if [ -z "$gt" ]; then
        printf "%-60s %s\n" "$name" "NO GT FOUND"
        continue
    fi

    read matches total <<< $(compare_accuracy "$pred" "$gt")
    printf "%-60s %s/%s\n" "$name" "$matches" "$total"
done

# --- Large-circuit bench dirs ---
echo
echo "===== Large-Circuit Confirmation (bench dirs) ====="
for bench_dir in "$PROJECT_DIR"/bench_obs_better_binding/bench_* "$PROJECT_DIR"/bench_36obs_d21_p001_2058r_100s; do
    if [ -d "$bench_dir/preds" ]; then
        echo
        echo "--- $(basename "$bench_dir") ---"
        "$PROJECT_DIR/scripts/compare_preds.sh" "$bench_dir"
    fi
done
