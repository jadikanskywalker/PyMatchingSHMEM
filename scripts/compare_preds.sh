#!/bin/bash
# Compare each final prediction file in a bench preds/ dir against ground truth.
# Usage: ./compare_preds.sh <bench_dir>
# Example: ./compare_preds.sh bench_36obs_d21_p001_2058r_100s

set -euo pipefail

BENCH_DIR="${1:?Usage: $0 <bench_dir>}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

BENCH_DIR="$(cd "$BENCH_DIR" 2>/dev/null && pwd || echo "$PROJECT_DIR/$BENCH_DIR")"
PREDS_DIR="$BENCH_DIR/preds"

if [ ! -d "$PREDS_DIR" ]; then
    echo "Error: $PREDS_DIR not found" >&2
    exit 1
fi

# Derive ground truth path from bench dir name:
#   bench_36obs_d21_p001_2058r_100s -> actual_obs_flips_36obs_d21_p001_2058r_100s.01
BENCH_NAME="$(basename "$BENCH_DIR")"
SUFFIX="${BENCH_NAME#bench_}"
GT_FILE="$PROJECT_DIR/testdems/actual_obs_flips_${SUFFIX}.01"

if [ ! -f "$GT_FILE" ]; then
    echo "Error: ground truth not found at $GT_FILE" >&2
    exit 1
fi

GT_LINES=$(wc -l < "$GT_FILE")
GT_BITS=$(awk '{ n += length($0) } END { print n }' "$GT_FILE")

printf "%-65s %s/%s\n" "FILE" "MATCH" "TOTAL"
printf "%0.s-" {1..80}; echo

for pred in "$PREDS_DIR"/preds_*.01; do
    name="$(basename "$pred")"

    # Skip per-PE files (intermediate outputs)
    if [[ "$name" == *_pe[0-9]* ]]; then
        continue
    fi

    match=$(paste -d "" "$pred" "$GT_FILE" | awk -v bpl="${GT_BITS}" '
    {
        n = length($0) / 2
        for (i = 1; i <= n; i++) {
            if (substr($0, i, 1) == substr($0, n + i, 1)) c++
        }
    }
    END { print c + 0 }')

    printf "%-65s %s/%s\n" "$name" "$match" "$GT_BITS"
done
