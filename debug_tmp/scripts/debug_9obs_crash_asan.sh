#!/bin/bash
#SBATCH --job-name=debug_9obs_asan
#SBATCH --output=debug_tmp/out/debug_9obs_asan-%j.out
#SBATCH --error=debug_tmp/out/debug_9obs_asan-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=128GB

# Small/cheap complement to debug_36obs_crash_gdb.sh: 9obs DEM at 2 PEs. Star topology (hub=obs8)
# splits 9 obs as PE0=[0-4] PE1=[5-8] (base/rem division), so obs8 (the hub, on PE1) connects to
# both local seams (obs5-7, also PE1) and CRTs (obs0-4, PE0<->PE1) -- same local-seam+CRT mix as
# the 36obs double-free, at a fraction of the size/time, under build_sos_debug (ASan+UBSan+SHMEM).

source ~/.bash_profile
conda activate pymatching

M=21
k=1
dem=~/PyMatchingSHMEM/testdems/error_model_9obs_d7_p01_337r.dem
det=~/PyMatchingSHMEM/testdems/detection_events_9obs_d7_p01_337r_50s.b8

outdir=~/PyMatchingSHMEM/debug_tmp/out/debug_9obs_asan_$SLURM_JOB_ID
mkdir -p "$outdir"
cd "$outdir"

export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=4G
export OMP_NUM_THREADS=16
export OMP_PLACES="cores(16)"
export OMP_PROC_BIND=true
export ASAN_OPTIONS=abort_on_error=1:halt_on_error=1

oshrun \
    -n 2 \
    --map-by ppr:1:package:PE=16 \
    --bind-to core \
    --report-bindings \
    ~/PyMatchingSHMEM/scripts/pe_output_wrapper.sh \
    ~/PyMatchingSHMEM/build_sos_debug/pymatching predict \
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
