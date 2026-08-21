#!/bin/bash
#SBATCH --job-name=debug_36obs_prte_st
#SBATCH --output=debug_tmp/out/debug_36obs_prte_st-%j.out
#SBATCH --error=debug_tmp/out/debug_36obs_prte_st-%j.err
#SBATCH --partition=zen4
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=768GB

# Reproduces the ntasks=2/sockets=2/ntps=1/nthreads=64 crash from the 36obs benchmark sweep
# (originally job 144389: rank 0 exited on signal 6 == the glibc "double free or corruption"
# SIGABRT), this time using PRTE's own built-in crash stack-trace facility instead of gdb.
# PRTE intercepts signals 6,7,8,11 (SIGABRT/SIGBUS/SIGFPE/SIGSEGV) by default
# (MCA prte_signal="6,7,8,11") and can dump a per-rank stack trace before the job dies -- no
# ptrace attach race, no zombie-process window, no sustained debugger overhead to fight PRTE's
# own liveness timeout with (see gdb-based attempts: full-lifetime wrap got killed as "hung",
# late-attach lost the race to the zombie-reap window). This is the native mechanism gdb was a
# substitute for; it was never explicitly enabled/pointed at a file in prior runs.

source ~/.bash_profile
conda activate pymatching

M=21
k=1
dem=~/PyMatchingSHMEM/testdems/error_model_36obs_d21_p001_2058r.dem
det=~/PyMatchingSHMEM/testdems/detection_events_36obs_d21_p001_2058r_100s.b8

outdir=~/PyMatchingSHMEM/debug_tmp/out/debug_36obs_prte_st_$SLURM_JOB_ID
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
    --mca prte_signal 6,7,8,11 \
    --mca prte_stacktrace_output file:stacktrace \
    --mca prte_timeout_for_stack_trace 60 \
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
