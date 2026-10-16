#!/bin/bash
#SBATCH --job-name=gen_128obs_m22
#SBATCH --output=dem_gen/out/gen_128obs_m22-%j.out
#SBATCH --partition=zen4
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=1400GB

# M=22 counterpart to gen_128obs.sh: preset_128obs()'s seam schedule was generated with its
# own internal M=22 (namespaced_main.cc calls preset_128obs() with no args, so it uses the
# default), but gen_128obs.sh decoded it at --rounds_per_partition 21, offsetting every seam
# from its containing partition's boundary. This version decodes at --rounds_per_partition 22
# instead, so partitions land exactly on the seam schedule's own boundaries -- both are kept
# (distinct round counts -> distinct cache filenames, no collision) to compare whether
# misaligned seams cost any decode performance (or, per the plan this session, to have a
# fallback in case the task-tree generality work doesn't yet handle mid-partition seams
# cleanly under M=21).
#
# --gen_rounds 1474: preset_128obs()'s own exact schedule length is 1472 (see gen_128obs.sh),
# rounded up to the next multiple of 22 (=22*67) instead of 21 this time.
#
# See gen_72obs.sh for the -n 1 (not 2) rationale.
#
# --use_threads IS required here (see gen_72obs.sh) -- it gates whether
# detector_error_model_to_user_graph actually populates virtual_boundaries/num_partitions,
# which is exactly what gets serialized into the cache. Omitting it silently caches a
# zero-partition graph. --mem bumped to 1400GB to cover the added threaded-decode sanity
# check this now runs.

cd ~/PyMatchingSHMEM/testdems
source ~/.bash_profile
conda activate pymatching
export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=128G

d=21
p=0.001
rounds=1474
shots=100
p_tag=$(echo "$p" | sed 's/^0\.//')
config_name="128obs_d${d}_p${p_tag}_${rounds}r"
graph_cache_name=graph_${config_name}_M22.cache
det_name=detection_events_${config_name}_${shots}s.b8
flips_name=actual_obs_flips_${config_name}_${shots}s.01
echo $graph_cache_name
echo $det_name
echo $flips_name

oshrun -n 1 \
    ~/PyMatchingSHMEM/build_sos/pymatching predict \
    --gen_code surface_code \
    --gen_task rotated_memory_x \
    --gen_distance $d \
    --gen_rounds $rounds \
    --gen_depolarization $p \
    --gen_surgery_preset 128obs \
    --gen_num_obs 128 \
    --graph_cache_path $graph_cache_name \
    --gen_det_out $det_name \
    --gen_obs_out $flips_name \
    --gen_sample_shots $shots \
    --in $det_name --in_format b8 \
    --out out_gen_sanity.01 --out_format 01 \
    --rounds_per_partition 22 \
    --task_division_strategy observable \
    --use_threads

echo "Done"
