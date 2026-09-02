#!/bin/bash
#SBATCH --job-name=scorep_variance_sweep
#SBATCH --output=debug_tmp/out/scorep_variance_sweep-%j.out
#SBATCH --error=debug_tmp/out/scorep_variance_sweep-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G

# Leaf/fusion decode-time variance (Score-P OTF2 trace) + sibling-steal count characterization
# (Score-P profile) at p=0.001 vs p=0.01. Same DEM/M/distance conventions as
# scripts/profile_threads_sweep.sh (d=21, surface_code, rotated_memory_x, rounds=10000, M=32) for
# direct comparability with the existing run_d21_profile_* data.
#
# Two shot counts per p: 100 shots for the trace job (keeps OTF2 trace size proportionate) and
# 1000 shots for the profile-only sibling-steal sweep (cheap to collect, better statistics for a
# rare event).

source ~/.bash_profile
conda activate pymatching

dir=~/PyMatchingSHMEM/debug_tmp/out/scorep_variance
mkdir -p $dir
cd $dir

d=21
rounds=10000
M=32
p_decs=(001 01)
thread_sweep=(8 16 32 64)

for p_dec in "${p_decs[@]}"; do
    mkdir -p p${p_dec}
    cd p${p_dec}

    stim gen \
        --rounds=$((rounds-1)) \
        --distance=$d \
        --after_clifford_depolarization=0.${p_dec} \
        --code surface_code \
        --task rotated_memory_x \
        > circuit.stim
    stim analyze_errors \
        --decompose_errors \
        --fold_loops \
        --in circuit.stim \
        > error_model.dem

    stim detect \
        --in circuit.stim \
        --shots 100 \
        --obs_out actual_obs_flips_100.01 \
        --obs_out_format 01 \
        --out detection_events_100.b8 \
        --out_format b8
    stim detect \
        --in circuit.stim \
        --shots 1000 \
        --obs_out actual_obs_flips_1000.01 \
        --obs_out_format 01 \
        --out detection_events_1000.b8 \
        --out_format b8

    # Job 1: leaf/fusion decode-time variance (trace mode), one representative thread count --
    # extraction_slot_wait was already substantial at 32 threads in the existing data.
    mkdir -p log_trace_32threads
    sbatch \
        --chdir=log_trace_32threads \
        --output=../log_trace_32threads.out \
        ~/PyMatchingSHMEM/debug_tmp/scripts/scorep_variance_call.sh \
            trace 32 $(pwd)/detection_events_100.b8 $(pwd)/error_model.dem

    # Job 2: sibling-steal counts (profile mode only), full thread sweep.
    for t in "${thread_sweep[@]}"; do
        mkdir -p log_profile_${t}threads
        sbatch \
            --chdir=log_profile_${t}threads \
            --output=../log_profile_${t}threads.out \
            ~/PyMatchingSHMEM/debug_tmp/scripts/scorep_variance_call.sh \
                profile $t $(pwd)/detection_events_1000.b8 $(pwd)/error_model.dem
    done

    cd ..
done
