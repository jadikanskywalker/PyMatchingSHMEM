#!/bin/bash
#SBATCH --job-name=test_65obs_with_seams_shmem_2pe
#SBATCH --output=debug_tmp/out/sbatch/test_65obs_with_seams_shmem_2pe-%j.out
#SBATCH --error=debug_tmp/out/sbatch/test_65obs_with_seams_shmem_2pe-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=64GB

# SHMEM/2-PE counterpart to test_65obs_with_seams.sh -- same 65-observable DEM and 16-gate
# --gen_surgery_spec (so some of those gates now land cross-PE -> CrossRankTask, not just
# LocalSeamTask), decoded across 2 PEs via build_sos.
#
# build_sos is compiled with USE_SHMEM, so the binary calls shmem_init_thread() unconditionally
# for every `predict` invocation (namespaced_main.cc's pm::main, not gated on --use_threads) --
# even pure DEM generation and the serial baseline must run under oshrun (-n 1), or the process
# aborts with no PMIx context. Also: under USE_SHMEM, --out always gets "_pe<pid>" inserted before
# the extension (namespaced_main.cc:101-121), unconditionally -- even at -n 1 -- so every --out
# path below is read back with that suffix.

source ~/.bash_profile
conda activate pymatching

cd ~/PyMatchingSHMEM

outdir=~/PyMatchingSHMEM/debug_tmp/out/test_65obs_with_seams_shmem_2pe_$SLURM_JOB_ID
mkdir -p "$outdir"
cd "$outdir"

BIN=~/PyMatchingSHMEM/build_sos/pymatching

export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=16G

# --- Generate the DEM + sample shots (-n 1: --gen_det_out/--gen_obs_out/--dem_cache_path are
# plain paths, not PE-suffixed, but the binary itself still needs a PMIx context to start) ---
export OMP_NUM_THREADS=1
oshrun -n 1 "$BIN" predict \
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
    --out gen_discard.01 \
    --out_format 01 \
    > gen.out 2>&1
echo "gen exit=$?"
grep "Generated DEM" gen.out

# --- Serial baseline (-n 1, no --use_threads) -- --out gets "_pe0" inserted regardless ---
oshrun -n 1 "$BIN" predict \
    --dem model.dem \
    --in det.b8 \
    --in_format b8 \
    --out predicted_serial.01 \
    --out_format 01 \
    > log_serial.out 2>log_serial.err
echo "serial exit=$?"
cp predicted_serial_pe0.01 predicted_serial.01

# --- Threaded DecodingUnit path across 2 PEs. my_obs_count per PE is ceil(65/2)=33 max, so
# PE=40 threads/PE covers it with margin. Some of the 16 surgery-spec gates now straddle the
# PE boundary (obs-major contiguous split -> PE0 owns obs 0-32, PE1 owns obs 33-64) and become
# CrossRankTasks instead of LocalSeamTasks. ---
export OMP_NUM_THREADS=40
export OMP_PLACES=cores
export OMP_PROC_BIND=true
oshrun \
    -n 2 \
    --map-by ppr:1:package:PE=40 \
    --bind-to core \
    --report-bindings \
    ~/PyMatchingSHMEM/scripts/pe_output_wrapper.sh \
    "$BIN" predict \
        --dem model.dem \
        --in det.b8 \
        --in_format b8 \
        --out predicted_threaded.01 \
        --out_format 01 \
        --rounds_per_partition 8 \
        --task_division_strategy observable \
        --cross_rank_fusion_window_size 1 \
        --use_threads
echo "threaded exit=$?"

python3 ~/PyMatchingSHMEM/scripts/combine_results.py predicted_threaded.01 2

echo
echo "Diff (serial vs combined threaded), should be empty:"
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
