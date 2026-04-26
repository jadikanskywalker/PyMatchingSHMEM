#!/bin/bash
#SBATCH --job-name=pymatching
#SBATCH --output=bench_obs.out
#SBATCH --error=bench_obs.err
#SBATCH --partition=zen4
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB

if [ $# -le 6 ]
  then
    echo "Args: [surgery_preset] [d] [p_dec] [shots] [rounds] [M] [k]"
    exit 1
else
    surgery_preset=$1
    d=$2
    p_dec=$3
    shots=$4
    rounds=$5
    M=($6)
    k=$7
fi

source ~/.bashrc
conda activate pymatching

if [ ! -d "bench_obs_${p_dec}" ]
  then
    mkdir bench_obs_${p_dec}
fi

cd bench_obs_${p_dec}

dem_suffix=${surgery_preset}_d${d}_p${p_dec}_${rounds}r
det_suffix=${dem_suffix}_${shots}s
dem=../testdems/error_model_$dem_suffix.dem
det=../testdems/detection_events_$det_suffix.b8
flips=../testdems/actual_obs_flips_$det_suffix.01
echo $dem
echo $det
echo $flips

shmem_threads=(16 32 64 128)
shmem_n1_n=(1 1 1 1)   # nodes
shmem_n1_pps=(1 1 1 1) # processes-per-socket
shmem_n1_nspn=(1 1 1 1)  # num sockets per node
shmem_n2_n=(1 1 1 1)
shmem_n2_pps=(2 2 2 1)
shmem_n2_nspn=(1 1 1 2)
shmem_n4_n=(1 1 1 2)
shmem_n4_pps=(4 4 2 1)
shmem_n4_nspn=(1 1 2 2)
shmem_n8_n=(1 1 2 4)
shmem_n8_pps=(8 4 2 1)
shmem_n8_nspn=(1 2 2 2)

serial_build=~/PyMatchingSHMEM/build/pymatching
# threads_build=~/PyMatchingSHMEM/build_threads_release/pymatching
shmem_build=~/PyMatchingSHMEM/build_sos/pymatching
echo "serial_build:  $serial_build" >> bench.out
# echo "threads_build: $threads_build" >> bench.out
echo "shmem_build: $shmem_build" >> bench.out
echo "----------" >> bench.out
echo "shots: $shots    rounds: $rounds" >> bench.out

for ((m=0; m<${#M[@]}; m++ )); do
    thisM=${M[$m]}
    echo "----------" >> bench.out
    echo "M: $thisM" >> bench.out
    echo "  SHMEM:" >> bench.out
    for ((i=0; i<${#shmem_threads[@]}; i++ )); do
        thisThreads=${shmem_threads[$i]}

        # thisN=${shmem_n1_n[$i]}
        # thisPPN=${shmem_n1_ppn[$i]}
        # sbatch \
        #     --nodes=$thisN \
        #     --ntasks-per-node=$thisPPN \
        #     --cpus-per-task=$thisThreads \
        #     ../scripts/benchmark_shmem_obs_call.sh \
        #         $thisN $thisPPN $thisThreads $thisM $k $dem $det $flips

        thisN=${shmem_n2_n[$i]}
        thisPPS=${shmem_n2_pps[$i]}
        thisNSPN=${shmem_n2_nspn[$i]}
        thisPPN=$((thisPPS*thisNSPN))
        thisMEM=$((thisPPN * 64))
        sbatch \
            --nodes=$thisN \
            --ntasks-per-node=$thisPPN \
            --cpus-per-task=$thisThreads \
            --mem=${thisMEM}GB \
            ../scripts/benchmark_shmem_obs_call.sh \
                $thisN $thisPPS $thisNSPN $thisThreads $thisM $k $dem $det $flips

        thisN=${shmem_n4_n[$i]}
        thisPPS=${shmem_n4_pps[$i]}
        thisNSPN=${shmem_n4_nspn[$i]}
        thisPPN=$((thisPPS*thisNSPN))
        thisMEM=$((thisPPN * 64))
        sbatch \
            --nodes=$thisN \
            --ntasks-per-node=$thisPPN \
            --cpus-per-task=$thisThreads \
            --mem=${thisMEM}GB \
            ../scripts/benchmark_shmem_obs_call.sh \
                $thisN $thisPPS $thisNSPN $thisThreads $thisM $k $dem $det $flips

        thisN=${shmem_n8_n[$i]}
        thisPPS=${shmem_n8_pps[$i]}
        thisNSPN=${shmem_n8_nspn[$i]}
        thisPPN=$((thisPPS*thisNSPN))
        thisMEM=$((thisPPN * 64))
        sbatch \
            --nodes=$thisN \
            --ntasks-per-node=$thisPPN \
            --cpus-per-task=$thisThreads \
            --mem=${thisMEM}GB \
            ../scripts/benchmark_shmem_obs_call.sh \
                $thisN $thisPPS $thisNSPN $thisThreads $thisM $k $dem $det $flips
    done
done

# start_serial=$(date +%s)
# $serial_build predict \
#     --dem $dem \
#     --in $det\
#     --in_format b8 \
#     --out preds_0.01 \
#     --out_format 01 \
#     > log_0.out
# end_serial=$(date +%s)
# serial_time=$((end_serial - start_serial))
# echo "0: $serial_time seconds" >> bench.out