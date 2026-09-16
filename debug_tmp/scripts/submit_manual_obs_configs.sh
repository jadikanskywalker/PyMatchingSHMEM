#!/bin/bash
# Reusable helper functions for manually (re)submitting SPECIFIC many-observable SHMEM benchmark
# configs -- e.g. resubmitting just the configs that failed in a prior sweep, or a one-off gdb
# debug pass over a preset -- instead of resubmitting benchmark_shmem_obs(_cached).sh's whole
# sweep. See memory project_multiobs_benchmark_pipeline for the pipeline this plugs into, and
# project_shmem_crash_debug_gdb for the *_gdb call-script variants used below.
#
# This is a function library, not a standalone job: source it, then call submit_dem / submit_cached
# / submit_dem_gdb / submit_cached_gdb per config. See the worked examples at the bottom (comment
# out/edit before running -- this is meant to be adapted per debugging session, same as every other
# debug_tmp/scripts/*.sh).
#
# All submissions assume dual-socket, 128-core/socket nodes and --exclusive (so --mem=1490GB is
# always safe/cheap regardless of preset size -- see benchmark_shmem_obs_cached.sh's own comment).
# `cd`s into the preset's own bench_<det_suffix>/ dir (created if missing) before each sbatch so
# the call script's relative ../preds, ../bench.out, out/ paths land in the same place the normal
# driver scripts use.

NODE_MEM=1490GB

# submit_dem preset rounds M nodes sockets ntpn ntasks nthreads [parent_dir]
submit_dem() {
    local preset=$1 rounds=$2 M=$3 nodes=$4 sockets=$5 ntpn=$6 ntasks=$7 nthreads=$8
    local parent_dir=${9:-.}
    local dem_suffix=${preset}_d21_p001_${rounds}r
    local det_suffix=${dem_suffix}_100s
    local d="$parent_dir/bench_$det_suffix"
    mkdir -p "$d/preds" "$d/out"
    local dem=~/PyMatchingSHMEM/testdems/error_model_$dem_suffix.dem
    local det=~/PyMatchingSHMEM/testdems/detection_events_$det_suffix.b8
    local flips=~/PyMatchingSHMEM/testdems/actual_obs_flips_$det_suffix.01
    ( cd "$d" && sbatch \
        --nodes=$nodes --exclusive --mem=$NODE_MEM \
        ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call.sh \
            $ntasks $sockets $ntpn $nthreads $M 1 $dem $det $flips )
}

# submit_cached preset rounds M nodes sockets ntpn ntasks nthreads cache_suffix [parent_dir]
# cache_suffix: "" for the plain M21 caches (72/144obs), "_M22" for the M22 caches (64/128/256obs)
submit_cached() {
    local preset=$1 rounds=$2 M=$3 nodes=$4 sockets=$5 ntpn=$6 ntasks=$7 nthreads=$8 suffix=$9
    local parent_dir=${10:-.}
    local dem_suffix=${preset}_d21_p001_${rounds}r
    local det_suffix=${dem_suffix}_100s
    local d="$parent_dir/bench_$det_suffix"
    mkdir -p "$d/preds" "$d/out"
    local graph_cache=~/PyMatchingSHMEM/testdems/graph_${dem_suffix}${suffix}.cache
    local det=~/PyMatchingSHMEM/testdems/detection_events_$det_suffix.b8
    local flips=~/PyMatchingSHMEM/testdems/actual_obs_flips_$det_suffix.01
    ( cd "$d" && sbatch \
        --nodes=$nodes --exclusive --mem=$NODE_MEM \
        ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call_cached.sh \
            $ntasks $sockets $ntpn $nthreads $M 1 $graph_cache $det $flips )
}

# Same args as submit_dem/submit_cached, but through the gdb-wrapped call scripts (pe_gdb_wrapper.sh
# instead of pe_output_wrapper.sh) -- catches a full "thread apply all bt full" backtrace on
# SIGSEGV/SIGABRT. See project_shmem_crash_debug_gdb for the known PRTE-liveness-kill risk.
submit_dem_gdb() {
    local preset=$1 rounds=$2 M=$3 nodes=$4 sockets=$5 ntpn=$6 ntasks=$7 nthreads=$8
    local parent_dir=${9:-.}
    local dem_suffix=${preset}_d21_p001_${rounds}r
    local det_suffix=${dem_suffix}_100s
    local d="$parent_dir/bench_$det_suffix"
    mkdir -p "$d/preds" "$d/out"
    local dem=~/PyMatchingSHMEM/testdems/error_model_$dem_suffix.dem
    local det=~/PyMatchingSHMEM/testdems/detection_events_$det_suffix.b8
    local flips=~/PyMatchingSHMEM/testdems/actual_obs_flips_$det_suffix.01
    ( cd "$d" && sbatch \
        --nodes=$nodes --exclusive --mem=$NODE_MEM \
        ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call_gdb.sh \
            $ntasks $sockets $ntpn $nthreads $M 1 $dem $det $flips )
}

submit_cached_gdb() {
    local preset=$1 rounds=$2 M=$3 nodes=$4 sockets=$5 ntpn=$6 ntasks=$7 nthreads=$8 suffix=$9
    local parent_dir=${10:-.}
    local dem_suffix=${preset}_d21_p001_${rounds}r
    local det_suffix=${dem_suffix}_100s
    local d="$parent_dir/bench_$det_suffix"
    mkdir -p "$d/preds" "$d/out"
    local graph_cache=~/PyMatchingSHMEM/testdems/graph_${dem_suffix}${suffix}.cache
    local det=~/PyMatchingSHMEM/testdems/detection_events_$det_suffix.b8
    local flips=~/PyMatchingSHMEM/testdems/actual_obs_flips_$det_suffix.01
    ( cd "$d" && sbatch \
        --nodes=$nodes --exclusive --mem=$NODE_MEM \
        ~/PyMatchingSHMEM/scripts/benchmark_shmem_obs_call_cached_gdb.sh \
            $ntasks $sockets $ntpn $nthreads $M 1 $graph_cache $det $flips )
}

# --- Worked examples (edit/uncomment as needed; this file is not meant to be run as-is) ---
#
# cd ~/PyMatchingSHMEM/bench_shmem_obs_2026_09_10   # match parent_dir to an existing sweep's output
#
# # Resubmit one specific failed config (e.g. after a placement-flag fix):
# submit_dem 36obs 2058 21 2 4 1 2 256
#
# # Full node/sockets/ntpn table per thread count (dual-socket 128-core/socket nodes), used to
# # resubmit/replay a whole preset's ntasks{1,2,4,8} x threads{16,32,64,128,256} sweep:
# threads=(16 32 64 128 256)
# n2_nodes=(1 1 1 1 2);   n2_sockets=(1 1 1 2 4);   n2_ntpn=(2 2 2 2 1)
# n4_nodes=(1 1 1 2 4);   n4_sockets=(1 1 2 4 8);   n4_ntpn=(4 4 4 2 1)
# n8_nodes=(1 1 2 4 8);   n8_sockets=(1 2 4 8 16);  n8_ntpn=(8 8 4 2 1)
# for i in "${!threads[@]}"; do
#     t=${threads[$i]}
#     submit_cached_gdb 72obs 2751 21 1 1 1 1 $t ""
#     submit_cached_gdb 72obs 2751 21 ${n2_nodes[$i]} ${n2_sockets[$i]} ${n2_ntpn[$i]} 2 $t ""
#     submit_cached_gdb 72obs 2751 21 ${n4_nodes[$i]} ${n4_sockets[$i]} ${n4_ntpn[$i]} 4 $t ""
#     submit_cached_gdb 72obs 2751 21 ${n8_nodes[$i]} ${n8_sockets[$i]} ${n8_ntpn[$i]} 8 $t ""
# done
