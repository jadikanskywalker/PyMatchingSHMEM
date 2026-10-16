#!/bin/bash
#SBATCH --job-name=test_65obs_with_seams
#SBATCH --output=debug_tmp/out/sbatch/test_65obs_with_seams-%j.out
#SBATCH --error=debug_tmp/out/sbatch/test_65obs_with_seams-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=64GB

# Same 65-observable, no-surgery-preset base as test_65obs_no_seams.sh, but with 16
# --gen_surgery_spec gates introducing real cross-observable seams -- spread across the full
# observable index range (0-61) so the same spec is reusable later for a multi-PE build_sos run
# (some gates should land same-PE -> LocalSeamTask, some cross-PE -> CrossRankTask, depending on
# ntasks), and staggered across the 255-round window so they don't all fall in the same partition.
# For now: single PE, build_threads only, exercising LocalSeamTask/divide_vb across observables
# (not just within one observable's own K_p partitions).

source ~/.bash_profile
conda activate pymatching

cd ~/PyMatchingSHMEM

outdir=~/PyMatchingSHMEM/debug_tmp/out/test_65obs_with_seams_$SLURM_JOB_ID
mkdir -p "$outdir"
cd "$outdir"

BIN=~/PyMatchingSHMEM/build_threads/pymatching

# --- Generate the DEM + sample shots, with 16 cross-observable gates (seams), staggered across
# both observable index (low..high, for later multi-PE same-PE/cross-PE coverage) and round (so
# they don't all land in the same partition) ---
$BIN predict \
    --gen_code repetition_code \
    --gen_task memory \
    --gen_distance 5 \
    --gen_rounds 255 \
    --gen_num_obs 65 \
    --gen_depolarization 0.01 \
    --gen_surgery_spec "0,1,20,3;4,5,35,3;8,9,50,3;12,13,65,3;16,17,80,3;20,21,95,3;24,25,110,3;28,29,125,3;32,33,140,3;36,37,155,3;40,41,170,3;44,45,185,3;48,49,200,3;52,53,215,3;56,57,230,3;60,61,245,3" \
    --gen_det_out det.b8 \
    --gen_obs_out actual_obs.01 \
    --gen_sample_shots 50 \
    --gen_sample_seed 42 \
    --dem_cache_path model.dem \
    --in det.b8 \
    --in_format b8 \
    --out /dev/null \
    --out_format 01 \
    > gen.out 2>&1
echo "gen exit=$?"
grep "Generated DEM" gen.out

# --- Serial baseline (no --use_threads) ---
export OMP_NUM_THREADS=1
$BIN predict \
    --dem model.dem \
    --in det.b8 \
    --in_format b8 \
    --out predicted_serial.01 \
    --out_format 01 \
    > log_serial.out 2>log_serial.err
echo "serial exit=$?"

# --- Threaded DecodingUnit path (>64 obs => needs_search_flooder true, per-obs solver needs
# max_threads >= my_obs_count=65 on a single PE) ---
export OMP_NUM_THREADS=72
export OMP_PLACES=cores
export OMP_PROC_BIND=true
$BIN predict \
    --dem model.dem \
    --in det.b8 \
    --in_format b8 \
    --out predicted_threaded.01 \
    --out_format 01 \
    --rounds_per_partition 8 \
    --task_division_strategy observable \
    --use_threads \
    > log_threaded.out 2>log_threaded.err
echo "threaded exit=$?"

echo
echo "Diff (serial vs threaded), should be empty:"
diff predicted_serial.01 predicted_threaded.01 && echo "MATCH"

# 65-bit lines (one line per shot, all observables concatenated) -- compare bit-by-bit, not line-by-line.
total_bits=$(awk '{ n += length($0) } END { print n }' actual_obs.01)
echo
echo "Serial vs actual_obs.01 (bits correct / total):"
paste -d "" predicted_serial.01 actual_obs.01 | awk -v total="$total_bits" '
{
    n = length($0) / 2
    for (i = 1; i <= n; i++) {
        if (substr($0, i, 1) == substr($0, n + i, 1)) c++
    }
}
END { print c + 0 "/" total }'
echo
echo "Threaded vs actual_obs.01 (bits correct / total):"
paste -d "" predicted_threaded.01 actual_obs.01 | awk -v total="$total_bits" '
{
    n = length($0) / 2
    for (i = 1; i <= n; i++) {
        if (substr($0, i, 1) == substr($0, n + i, 1)) c++
    }
}
END { print c + 0 "/" total }'
