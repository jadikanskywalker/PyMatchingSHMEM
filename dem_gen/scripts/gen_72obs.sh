#!/bin/bash
#SBATCH --job-name=gen_72obs
#SBATCH --output=dem_gen/out/gen_72obs-%j.out
#SBATCH --partition=zen4
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=1400GB

# Same pattern as gen_36obsx4.sh: --graph_cache_path caches the compact binary UserGraph
# (graph_cache.h/.cc) instead of the raw text DEM. --gen_rounds 2751 is preset_72obs()'s
# own exact schedule length (computed directly from its SurgerySpec list: max(start_round+
# duration) across every gate it emits, default M=21) -- using fewer rounds would silently
# clamp/collapse its later gates onto the last available round instead of spacing them out
# as intended (this is exactly what happened to the existing 64obs DEM, generated with
# --gen_rounds 704 against a preset needing 1407).
#
# -n 1, not 2: an earlier attempt at -n 2 OOM'd at 512GB -- turns out DEM generation isn't
# rank-guarded, so BOTH ranks independently regenerated the full ~88.5M-detector DEM
# concurrently (confirmed: "Generated DEM: ..." printed twice, interleaved), roughly
# doubling peak memory for no benefit at this stage.
#
# --use_threads IS required here, even though this is "just" a cache-build run: it gates
# detector_error_model_to_user_graph's call to partition_nodes_by_obs_patch
# (user_graph.cc:639-645), which populates virtual_boundaries/num_partitions/node_part_id --
# exactly the fields write_user_graph_cache serializes. Omitting it (as an earlier version of
# this script did) silently writes a cache with zero partitions; --graph_cache_path readers
# then throw "SharedMatchingGraph requires at least one partition" instead of falling back to
# a raw DEM. --mem bumped to 1400GB (of ~1.5TB/node) to cover the added DecodingUnit/
# partition/region-arena setup this now exercises during the tacked-on sanity decode.

cd ~/PyMatchingSHMEM/testdems
source ~/.bash_profile
conda activate pymatching
export FI_VERBS_DEVICE_NAME="mlx5_2"
export SHMEM_SYMMETRIC_SIZE=128G

d=21
p=0.001
rounds=2751
shots=100
p_tag=$(echo "$p" | sed 's/^0\.//')
config_name="72obs_d${d}_p${p_tag}_${rounds}r"
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
    --gen_surgery_preset 72obs \
    --gen_num_obs 73 \
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
