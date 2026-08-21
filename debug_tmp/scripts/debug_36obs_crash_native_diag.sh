#!/bin/bash
#SBATCH --job-name=debug_36obs_native_diag
#SBATCH --output=debug_tmp/out/debug_36obs_native_diag-%j.out
#SBATCH --error=debug_tmp/out/debug_36obs_native_diag-%j.err
#SBATCH --partition=zen4
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=768GB

# Reproduces the ntasks=2/sockets=2/ntps=1/nthreads=64 double-free crash from the 36obs benchmark
# sweep (originally job 144389: rank 0 exited on signal 6). Plain native run -- no gdb (attach
# races/zombie-reap window; full-lifetime wrap gets killed by PRTE's own liveness timeout under
# ptrace overhead), no PRTE --mca stacktrace flags (confirmed not to fire for SOS/OpenSHMEM
# binaries, which don't hook into OPAL's MPI_INIT-time signal handler install). Instead this
# relies on the new double-del diagnostic added directly to SHMEMArena<T>::del()
# (shmem_arena.h) -- GraphFillRegion's own arena, which is what job 144398's full gdb dump showed
# actually double-freed inside Mwpm::shatter_blossom_and_extract_matches.

source ~/.bash_profile
conda activate pymatching

M=21
k=1
dem=~/PyMatchingSHMEM/testdems/error_model_36obs_d21_p001_2058r.dem
det=~/PyMatchingSHMEM/testdems/detection_events_36obs_d21_p001_2058r_100s.b8

outdir=~/PyMatchingSHMEM/debug_tmp/out/debug_36obs_native_diag_$SLURM_JOB_ID
mkdir -p "$outdir"
cd "$outdir"

export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=16G
export OMP_NUM_THREADS=64
export OMP_PLACES="cores(64)"
export OMP_PROC_BIND=true

oshrun \
    -n 2 \
    --map-by ppr:1:package:PE=64 \
    --bind-to core \
    --report-bindings \
    ~/PyMatchingSHMEM/scripts/pe_output_wrapper.sh \
    ~/PyMatchingSHMEM/build_sos/pymatching predict \
        --dem "$dem" \
        --in "$det" \
        --in_format b8 \
        --out predicted.01 \
        --out_format 01 \
        --rounds_per_partition $M \
        --cross_rank_fusion_window_size $k \
        --task_division_strategy observable \
        --extraction_unit_size 4 \
        --extract_preemptively \
        --use_threads
