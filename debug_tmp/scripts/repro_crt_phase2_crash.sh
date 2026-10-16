#!/bin/bash
#SBATCH --job-name=repro_crt_phase2_crash
#SBATCH --output=debug_tmp/out/sbatch/repro_crt_phase2_crash-%j.out
#SBATCH --error=debug_tmp/out/sbatch/repro_crt_phase2_crash-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=2
#SBATCH --exclusive
#SBATCH --mem=1490GB
#SBATCH --exclude=rpc-96-1

# Reproduces the free(): invalid pointer / SIGABRT seen on rank 3 of the 36obs 4-task CRT
# validation run (crt_validation_runs/segfault_gdb_repro*), deterministic at shot 9, using
# pe_gdb_wrapper_crashing_only.sh to capture just the crashing thread's own backtrace first
# (the full "thread apply all" dump truncated at 60/146 threads before reaching it last time).

source ~/.bash_profile
conda activate pymatching

outdir=~/PyMatchingSHMEM/debug_tmp/out/repro_crt_phase2_crash_$SLURM_JOB_ID
mkdir -p "$outdir"
cd "$outdir"

export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=256G
export OMP_NUM_THREADS=128
export OMP_PLACES="cores(128)"
export OMP_PROC_BIND=true

oshrun \
    -n 4 \
    --map-by ppr:2:node:PE=128 \
    --bind-to core \
    --report-bindings \
    ~/PyMatchingSHMEM/debug_tmp/scripts/pe_gdb_wrapper_crashing_only.sh \
    ~/PyMatchingSHMEM/build_sos/pymatching predict \
        --graph_cache_path "$WORK/testdems/graph_36obs_d21_p001_2058r.cache" \
        --in "$WORK/testdems/detection_events_36obs_d21_p001_2058r_100s.b8" \
        --in_format b8 \
        --out predicted.01 \
        --out_format 01 \
        --rounds_per_partition 21 \
        --cross_rank_fusion_window_size 1 \
        --task_division_strategy observable \
        --extraction_unit_size 16 \
        --extract_preemptively \
        --use_threads
