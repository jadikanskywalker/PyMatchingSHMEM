#!/bin/bash
#SBATCH --job-name=verify_dense_arena
#SBATCH --output=debug_tmp/out/sbatch/verify_dense_arena-%j.out
#SBATCH --error=debug_tmp/out/sbatch/verify_dense_arena-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=32G

# Step A of the DenseArena-port plan (plans/snappy-toasting-bee.md): confirm the dense-arena
# threads build (ENABLE_DENSE_REGION_ARENA=ON) produces byte-identical predictions to the baseline
# threads build (Arena, ENABLE_DENSE_REGION_ARENA=OFF) -- a pure allocator/layout swap must not
# change any decode result. Sweeps both ROUND and OBS division strategies and a couple of thread
# counts, since region allocation order interacts with partition/thread assignment.

source ~/.bash_profile
conda activate pymatching

baseline=~/PyMatchingSHMEM/build_threads/pymatching
dense=~/PyMatchingSHMEM/build_threads_dense/pymatching

outdir=~/PyMatchingSHMEM/debug_tmp/out/verify_dense_arena
mkdir -p "$outdir"
cd "$outdir"

export OMP_PLACES=cores
export OMP_PROC_BIND=close

fail=0

run_and_compare() {
    local tag=$1 strategy=$2 dem=$3 det=$4 nthreads=$5 rounds_per_partition=$6
    export OMP_NUM_THREADS=$nthreads
    local dir="${tag}_t${nthreads}"
    mkdir -p "$dir"

    "$baseline" predict \
        --dem "$dem" --in "$det" --in_format b8 \
        --out "${dir}/baseline.01" --out_format 01 \
        --task_division_strategy "$strategy" \
        --rounds_per_partition "$rounds_per_partition" \
        --use_threads \
        > "${dir}/baseline.out" 2> "${dir}/baseline.err"
    local rc1=$?

    "$dense" predict \
        --dem "$dem" --in "$det" --in_format b8 \
        --out "${dir}/dense.01" --out_format 01 \
        --task_division_strategy "$strategy" \
        --rounds_per_partition "$rounds_per_partition" \
        --use_threads \
        > "${dir}/dense.out" 2> "${dir}/dense.err"
    local rc2=$?

    if [ $rc1 -ne 0 ] || [ $rc2 -ne 0 ]; then
        echo "${dir}: FAIL (exit codes baseline=$rc1 dense=$rc2 -- see ${dir}/*.err)"
        fail=1
        return
    fi

    if cmp -s "${dir}/baseline.01" "${dir}/dense.01"; then
        echo "${dir}: MATCH"
    else
        echo "${dir}: MISMATCH -- predictions differ"
        fail=1
    fi
}

# --- ROUND strategy: single-observable DEM (ROUND partitioning requires this -- see
# feedback/project memory: multi-obs DEMs misread the observable coordinate as round coordinate) ---
round_dem=~/PyMatchingSHMEM/acc_results/out/dems/1obs_d9_p0005_r200.dem
round_det=~/PyMatchingSHMEM/acc_results/out/events/1obs_d9_p0005_r200_1000s.b8
for t in 4 16; do
    run_and_compare round_strategy round "$round_dem" "$round_det" "$t" 20
done

# --- OBS strategy: 9-observable DEM ---
obs_dem=~/PyMatchingSHMEM/testdems/error_model_9obs_d21_p001_672r.dem
obs_det=~/PyMatchingSHMEM/testdems/detection_events_9obs_d21_p001_672r_100s.b8
for t in 9 18; do
    run_and_compare obs_strategy observable "$obs_dem" "$obs_det" "$t" 32
done

echo "=== overall: $([ $fail -eq 0 ] && echo PASS || echo FAIL) ==="
exit $fail
