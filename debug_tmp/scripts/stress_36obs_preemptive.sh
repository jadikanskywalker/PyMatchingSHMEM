#!/bin/bash
#SBATCH --job-name=stress_36obs_preemptive
#SBATCH --output=debug_tmp/out/stress_36obs_preemptive-%j.out
#SBATCH --error=debug_tmp/out/stress_36obs_preemptive-%j.err
#SBATCH --partition=zen4
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=768GB

# 36obs stress test for the decentralized shot-sync refactor (extraction_done / ShotIOResource):
# 2 PEs, 32 threads/PE, preemptive extraction, num_repeats=2 to also exercise reset() at full scale.

source ~/.bash_profile
conda activate pymatching

M=21
k=1
dem=~/PyMatchingSHMEM/testdems/error_model_36obs_d21_p001_2058r.dem
det=~/PyMatchingSHMEM/testdems/detection_events_36obs_d21_p001_2058r_100s.b8

outdir=~/PyMatchingSHMEM/debug_tmp/out/stress_36obs_preemptive_$SLURM_JOB_ID
mkdir -p "$outdir"
cd "$outdir"

export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=16G
export OMP_NUM_THREADS=32
export OMP_PLACES="cores(32)"
export OMP_PROC_BIND=true

start=$(date +%s)
oshrun \
    -n 2 \
    --map-by ppr:1:package:PE=32 \
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
        --use_threads \
        --num_repeats 2
end=$(date +%s)
echo "Run completed in $((end - start)) seconds."

python3 ~/PyMatchingSHMEM/scripts/combine_results.py predicted.01 2
wc -l predicted.01
