#!/bin/bash
#SBATCH --job-name=gen_144obs
#SBATCH --output=dem_gen/out/gen_144obs-%j.out
#SBATCH --partition=zen4
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=1450GB

# See gen_72obs.sh for the --graph_cache_path/--gen_rounds rationale. --gen_rounds 3444 is
# preset_144obs()'s own exact schedule length (two independent preset_72obs() copies ending
# at round 2751, +2*M idle gap, +one reduction_adder top-level carry at obs 146 -- mirrors
# preset_72obs's own construction from preset_36obs, one level higher; note --gen_num_obs is
# 147, not 144, matching preset_72obs's own 73-not-72 convention of adding a real carry
# observable). The largest of this batch: ~281GB transient DEM, 24108 partitions (more than
# gen_36obsx4.sh's 14436) -- --mem/SHMEM_SYMMETRIC_SIZE sized up accordingly.

# -n 1, not 2: an earlier attempt at -n 2 OOM'd -- DEM generation isn't rank-guarded, so
# BOTH ranks independently regenerated the full DEM concurrently, roughly doubling peak
# memory for no benefit at this stage.
#
# --use_threads IS required here (see gen_72obs.sh) -- it gates whether
# detector_error_model_to_user_graph actually populates virtual_boundaries/num_partitions,
# which is exactly what gets serialized into the cache. Omitting it silently caches a
# zero-partition graph. --mem bumped to 1450GB (largest preset in this batch, of ~1.5TB/node)
# to cover the added threaded-decode sanity check this now runs.

cd ~/PyMatchingSHMEM/testdems
source ~/.bash_profile
conda activate pymatching
export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=256G

d=21
p=0.001
rounds=3444
shots=100
p_tag=$(echo "$p" | sed 's/^0\.//')
config_name="144obs_d${d}_p${p_tag}_${rounds}r"
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
    --gen_surgery_preset 144obs \
    --gen_num_obs 147 \
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
