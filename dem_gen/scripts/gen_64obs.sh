#!/bin/bash
#SBATCH --job-name=gen_64obs
#SBATCH --output=dem_gen/out/gen_64obs-%j.out
#SBATCH --partition=zen4
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=1400GB

# Builds a --graph_cache_path cache for the existing (already correctly-sized, non-truncated)
# testdems/error_model_64obs_d21_p001_1408r.dem, at --rounds_per_partition 22 -- matching
# preset_64obs()'s own internal gate-timing M=22, so decode partitions land exactly on the
# preset's own seam-schedule boundaries (unlike the 128obs/256obs M=21 caches, whose seams
# were generated at M=22 internally but decoded at rounds_per_partition=21, offsetting seams
# from partition boundaries). Reuses the existing DEM/100-shot sample directly (--dem, not
# --gen_code) -- no need to regenerate, 1408 is already the correct round count for M=22
# (1407 strict minimum, rounded up to the next multiple of 22).
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
export SHMEM_SYMMETRIC_SIZE=32G

graph_cache_name=graph_64obs_d21_p001_1408r_M22.cache
det_name=detection_events_64obs_d21_p001_1408r_100s.b8
flips_name=actual_obs_flips_64obs_d21_p001_1408r_100s.01
echo $graph_cache_name
echo $det_name
echo $flips_name

# -n 1: see gen_72obs.sh notes -- generation/caching isn't rank-guarded, so multiple ranks
# would redundantly duplicate peak memory for no benefit at this stage.
oshrun -n 1 \
    ~/PyMatchingSHMEM/build_sos/pymatching predict \
    --dem error_model_64obs_d21_p001_1408r.dem \
    --graph_cache_path $graph_cache_name \
    --in $det_name --in_format b8 \
    --out out_gen_sanity.01 --out_format 01 \
    --rounds_per_partition 22 \
    --task_division_strategy observable \
    --use_threads

echo "Done"
