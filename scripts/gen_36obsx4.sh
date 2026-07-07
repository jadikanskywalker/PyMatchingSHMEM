#!/bin/bash
#SBATCH --job-name=gen_36obsx4
#SBATCH --output=gen_36obsx4-%j.out
#SBATCH --partition=zen4
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=1024GB

if [ $# -le 3 ]
  then
    echo "Args: [d] [p (e.g. 0.001)] [rounds] [shots]"
    exit 1
else
    d=$1
    p=$2
    rounds=$3
    shots=$4
fi

cd ~/PyMatchingSHMEM/testdems
source ~/.bash_profile
conda activate pymatching
export FI_VERBS_DEVICE_NAME="mlx5_2"
# See scripts/gen_36obsx4.sh notes: for this 133M-detector graph (d=21, 8400 rounds,
# 36 observables, 14436 partitions), each PE's symmetric heap must hold buffers sized
# by the FULL graph (they do NOT shrink with more PEs -- MatchingGraph is the whole
# global structure on every rank, only work is partitioned, not the node/region arrays):
#   node_ephemeral_fields_ptr: NUM_BUFFERS_PER_UNIT(2) * 133,347,568 * sizeof(DetectorNodeEphemeralFields=64)   ~= 15.9 GiB
#   regions_ptr:               regions_nelems_per_solver(9280) * num_partitions(14436) * 2 * sizeof(GraphFillRegion=144) ~= 35.9 GiB
#   child_edges_ptr:           same nelems * num_partitions * 2 * sizeof(BlossomChild=40)                        ~= 10.0 GiB
#   (task_status_ptr/task_fusion_summary_ptr scale with num_cross_rank_fusions, not node count -- negligible here)
# Total ~62 GiB/PE measured against the earlier failed run (16G heap, ~18G overrun on
# regions_ptr after node_ephemeral_fields_ptr had already consumed ~15.9G of it).
# 128G/PE gives ~2x margin per PE; with 2 PEs on one exclusive node (1.5TB RAM) this is
# comfortably within the node's actual memory.
export SHMEM_SYMMETRIC_SIZE=128G

p_tag=$(echo "$p" | sed 's/^0\.//')
config_name="36obsx4_d${d}_p${p_tag}_${rounds}r"
# Deliberately NOT using --dem_cache_path here: for a graph this size the raw text DEM
# is enormous (the p=0.01 attempt produced a 189GB .dem file). --graph_cache_path
# caches the compact binary UserGraph instead (see graph_cache.h/.cc), skipping DEM
# text entirely on any future rerun. The DEM is still built once in-memory this run
# (unavoidable -- that's the actual circuit/detector synthesis), it's just never
# serialized to text.
graph_cache_name=graph_$config_name.cache
det_name=detection_events_${config_name}_${shots}s.b8
flips_name=actual_obs_flips_${config_name}_${shots}s.01
echo $graph_cache_name
echo $det_name
echo $flips_name

# 2 PEs (one per socket on this node) -- a single-PE run wouldn't exercise cross-rank
# fusion at all, which is the main thing worth stress-testing on a graph this size.
oshrun -n 2 --map-by ppr:2:package:PE=32 --bind-to core --report-bindings \
    ~/PyMatchingSHMEM/build_sos/pymatching predict \
    --gen_code surface_code \
    --gen_task rotated_memory_x \
    --gen_distance $d \
    --gen_rounds $rounds \
    --gen_depolarization $p \
    --gen_surgery_preset 36obsx4 \
    --gen_num_obs 36 \
    --graph_cache_path $graph_cache_name \
    --gen_det_out $det_name \
    --gen_obs_out $flips_name \
    --gen_sample_shots $shots \
    --in $det_name --in_format b8 \
    --out out_gen_sanity.01 --out_format 01 \
    --rounds_per_partition 21 --obs_coors_included \
    --cross_rank_fusion_window_size 1 \
    --task_division_strategy observable --use_threads

echo "Done"
