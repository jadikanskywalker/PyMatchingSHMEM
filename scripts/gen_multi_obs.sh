#!/bin/bash
#SBATCH --job-name=gen_multi_obs
#SBATCH --output=gen_multi_obs-%j.out
#SBATCH --error=gen_multi_obs-%j.err
#SBATCH --partition=zen4
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=512GB

if [ $# -le 7 ]
  then
    echo "Args: [code] [task] [d] [0.p] [rounds] [surgery_preset] [num_obs] [shots]"
    exit 1
else
    code=$1
    task=$2
    d=$3
    p_dec=$4
    rounds=$5
    surgery_preset=$6
    num_obs=$7
    shots=$8
fi

cd ~/PyMatchingSHMEM/testdems
source ~/.bashrc
conda activate pymatching

config_name="${surgery_preset}_d${d}_p${p_dec}_${rounds}r"
dem_name=error_model_$config_name.dem
det_name=detection_events_${config_name}_${shots}s.b8
flips_name=actual_obs_flips_${config_name}_${shots}s.01
echo $dem_name
echo $det_name
echo $flips_name

p=0.$p_dec
rounds=$((rounds-1))

echo
echo Generating DEM
python3 ../scripts/gen_multi_obs.py \
    --num_observables $num_obs \
    --rounds $rounds \
    --distance $d \
    --after_clifford_depolarization $p \
    --code $code \
    --task $task \
    --surgery_preset $surgery_preset \
    --surgery_duration $d \
    > $dem_name
echo Sampling DEM
# Sample detection events FROM THE DEM (not the circuit) so seam errors fire
stim sample_dem \
    --in $dem_name \
    --shots $shots \
    --out $det_name \
    --out_format b8 \
    --obs_out $flips_name \
    --obs_out_format 01
echo Done

# python3 ../scripts/gen_multi_obs.py \
#     --num_observables 48 \
#     --rounds 50 \
#     --distance 5 \
#     --after_clifford_depolarization 0.01 \
#     --code repetition_code \
#     --task memory \
#     --surgery_preset 48obs \
#     --surgery_duration 2 \
#     > error_model.dem
