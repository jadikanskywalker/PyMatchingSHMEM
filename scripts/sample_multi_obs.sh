#!/bin/bash
#SBATCH --job-name=PyMatchingSHMEM
#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --partition=zen4
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB

if [ $# -le 7 ]
  then
    echo "Args: [name_append] [shots]"
    exit 1
else
    config_name=$1
    shots=$2
fi

config_name="${$surgery_preset}_d${$d}_p${$p}"
det_name=detection_events_$config_name.b8
flips_name=actual_obs_flips_$config_name.01
echo $dem_name
echo $det_name
echo $flips_name

stim sample_dem \
    --in $dem_name \
    --shots $shots \
    --out $det_name \
    --out_format b8 \
    --obs_out $flips_name \
    --obs_out_format 01