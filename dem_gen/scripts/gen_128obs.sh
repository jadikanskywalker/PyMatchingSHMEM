#!/bin/bash
#SBATCH --job-name=gen_128obs
#SBATCH --output=dem_gen/out/gen_128obs-%j.out
#SBATCH --partition=zen4
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=1400GB

# See gen_72obs.sh for the --graph_cache_path/--gen_rounds rationale. preset_128obs()'s own
# exact schedule length is 1472 (two independent preset_64obs() copies ending at round 1407,
# +2*M idle gap, +one plain top-level join gate -- mirrors preset_64obs's own sparser,
# non-adder top-join style rather than preset_72obs/144obs's denser reduction_adder style,
# deliberately, for comparing seam density cost). --gen_rounds must be an exact multiple of
# --rounds_per_partition (see project_dem_presets memory: a non-multiple round count breaks
# partition_nodes_by_obs_patch's uniform-K_p assumption) -- 1472 isn't a multiple of 21, so
# this rounds up to 1491 (=21*71), adding a few harmless idle rounds after the last real gate.

# -n 1, not 2: an earlier attempt at -n 2 OOM'd -- DEM generation isn't rank-guarded, so
# BOTH ranks independently regenerated the full DEM concurrently, roughly doubling peak
# memory for no benefit at this stage.
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
rounds=1491
shots=100
p_tag=$(echo "$p" | sed 's/^0\.//')
config_name="128obs_d${d}_p${p_tag}_${rounds}r"
graph_cache_name=graph_$config_name.cache
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
    --rounds_per_partition 21 \
    --task_division_strategy observable \
    --use_threads

echo "Done"
