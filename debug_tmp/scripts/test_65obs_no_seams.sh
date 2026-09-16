#!/bin/bash
#SBATCH --job-name=test_65obs_no_seams
#SBATCH --output=debug_tmp/out/sbatch/test_65obs_no_seams-%j.out
#SBATCH --error=debug_tmp/out/sbatch/test_65obs_no_seams-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=64GB

# First real exercise of the new SearchFlooder/>64-obs path: 65 fully independent observables
# (no --gen_surgery_preset/--gen_surgery_spec, so zero cross-observable gates -- no seams between
# observables at all), single PE, OBS task_division_strategy so each observable still gets K_p>1
# local partitions (exercises divide_vb/LocalSeamTask/extraction intra-observable, just not
# CrossRankTask/cross-PE fusion). Compares the threaded DecodingUnit path against the serial
# to_mwpm() path on the same detection events.

source ~/.bash_profile
conda activate pymatching

cd ~/PyMatchingSHMEM

outdir=~/PyMatchingSHMEM/debug_tmp/out/test_65obs_no_seams_$SLURM_JOB_ID
mkdir -p "$outdir"
cd "$outdir"

BIN=~/PyMatchingSHMEM/build_threads/pymatching

# --- Generate the DEM + sample shots, one shot (no surgery preset/spec => no seams) ---
$BIN predict \
    --gen_code repetition_code \
    --gen_task memory \
    --gen_distance 5 \
    --gen_rounds 255 \
    --gen_num_obs 65 \
    --gen_depolarization 0.001 \
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
export OMP_NUM_THREADS=128
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
    --extraction_unit_size 4 \
    --extract_preemptively \
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
