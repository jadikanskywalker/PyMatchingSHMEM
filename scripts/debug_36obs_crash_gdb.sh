#!/bin/bash
#SBATCH --job-name=debug_36obs_gdb
#SBATCH --output=debug_36obs_gdb-%j.out
#SBATCH --error=debug_36obs_gdb-%j.err
#SBATCH --partition=zen4
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=768GB

# Reproduces the ntasks=2/sockets=2/ntps=1/nthreads=64 double-free crash from the 36obs benchmark
# sweep. Wrapping the whole process lifetime in gdb (even with no MALLOC_CHECK_) was itself slow
# enough that PRTE's own liveness detection killed rank 0 as "hung" well before reaching the crash
# window -- so pe_gdb_late_attach_wrapper.sh runs natively and only attaches gdb ~12 minutes in,
# shortly before the crash window, for a brief one-time ptrace-attach pause instead of sustained
# overhead across the whole run.

source ~/.bash_profile
conda activate pymatching

M=21
k=1
dem=testdems/error_model_36obs_d21_p001_2058r.dem
det=testdems/detection_events_36obs_d21_p001_2058r_100s.b8

outdir=debug_36obs_gdb_$SLURM_JOB_ID
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
    ~/PyMatchingSHMEM/scripts/pe_gdb_late_attach_wrapper.sh \
    ~/PyMatchingSHMEM/build_sos/pymatching predict \
        --dem ../$dem \
        --in ../$det \
        --in_format b8 \
        --out predicted.01 \
        --out_format 01 \
        --rounds_per_partition $M \
        --cross_rank_fusion_window_size $k \
        --task_division_strategy observable \
        --extraction_unit_size 4 \
        --extract_preemptively \
        --use_threads
