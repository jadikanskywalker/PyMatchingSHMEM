#!/bin/bash
#SBATCH --job-name=pybatch
#SBATCH --output=bench_obs-%j.out
#SBATCH --error=bench_obs-%j.err
#SBATCH --partition=zen4
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=256GB

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

source ~/.bash_profile
conda activate pymatching

dem_suffix=${surgery_preset}_d${d}_p${p_dec}_${rounds}r
det_suffix=${dem_suffix}_${shots}s
dirname=bench_$det_suffix

if [ ! -d "$dirname" ]
  then
    mkdir $dirname
fi

cd $dirname

if [ ! -d preds ]
  then
    mkdir preds
fi

if [ ! -d out ]
  then
    mkdir out
fi

dem=../testdems/error_model_$dem_suffix.dem
det=../testdems/detection_events_$det_suffix.b8
flips=../testdems/actual_obs_flips_$det_suffix.01
echo $dem
echo $det
echo $flips


shmem_threads=(8 16 32 64)

# Fill sockets
shmem_ntasks2_sockets=(1 1 1 1 2)
shmem_ntasks2_ntps=(   2 2 2 2 1) # ntasks per socket

shmem_ntasks4_nodes=(  1 1 1 1 2)
shmem_ntasks4_sockets=(1 1 1 2 4)
shmem_ntasks4_ntps=(   4 4 4 2 1)

shmem_ntasks8_nodes=(  1 1 1 2 4)
shmem_ntasks8_sockets=(1 1 2 4 8)
shmem_ntasks8_ntps=(   8 8 4 2 1)

shmem_ntasks16_nodes=(   1 1 2 4 8)
shmem_ntasks16_sockets=( 1 2 4 8 16)
shmem_ntasks16_ntps=(   16 8 4 2 1)


# # Quarter subscribe sockets sockets
# shmem_quarter_threads=(8 16 16 32 32 32)
# shmem_quarter_ntasks=( 8  8  4  8  4  2)
# shmem_quarter_ntps=(   4  2  2  1  1  1)
# shmem_quarter_sockets=(2  4  2  8  4  2)
# shmem_quarter_nodes=(  1  2  1  4  2  1)

shmem_quarter_threads=(16 32 32)
shmem_quarter_ntasks=( 4  4  2)
shmem_quarter_ntps=(   2  1  1)
shmem_quarter_sockets=(2  4  2)
shmem_quarter_nodes=(  1  2  1)


# # Half subscribe sockets sockets
# shmem_half_threads=(16 32 32 64 64 64)
# shmem_half_ntasks=(  8  8  4  8  4  2)
# shmem_half_ntps=(    4  2  2  1  1  1)
# shmem_half_sockets=( 2  4  2  8  4  2)
# shmem_half_nodes=(   1  2  1  4  2  1)

shmem_half_threads=(32 64 64)
shmem_half_ntasks=(  4  4  2)
shmem_half_ntps=(    2  1  1)
shmem_half_sockets=( 2  4  2)
shmem_half_nodes=(   1  2  1)


repeats=1

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

        # # single socket runs
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=1 \
            --exclusive \
            --mem=256GB \
            ../scripts/benchmark_shmem_obs_call.sh \
                1 1 1 $thisThreads $thisM $k $dem $det $flips
        done

        thisSockets=${shmem_ntasks2_sockets[$i]}
        thisNTPS=${shmem_ntasks2_ntps[$i]}
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=1 \
            --exclusive \
            --mem=512GB \
            ../scripts/benchmark_shmem_obs_call.sh \
                2 $thisSockets $thisNTPS $thisThreads $thisM $k $dem $det $flips
        done

        thisNodes=${shmem_ntasks4_nodes[$i]}
        thisSockets=${shmem_ntasks4_sockets[$i]}
        thisSPN=$((thisSockets / thisNodes))
        thisNTPS=${shmem_ntasks4_ntps[$i]}
        thisNTPN=$((thisSPN * thisNTPS)) # ntasks per node
        thisMEM=$((4 * 160 / $thisNodes))
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=$thisNodes \
            --exclusive \
            --mem=${thisMEM}GB \
            ../scripts/benchmark_shmem_obs_call.sh \
                4 $thisSockets $thisNTPS $thisThreads $thisM $k $dem $det $flips
        done

        # thisNodes=${shmem_ntasks8_nodes[$i]}
        # thisSockets=${shmem_ntasks8_sockets[$i]}
        # thisSPN=$((thisSockets / thisNodes))
        # thisNTPS=${shmem_ntasks8_ntps[$i]}
        # thisNTPN=$((thisSPN * thisNTPS)) # ntasks per node
        # thisMEM=$((8 * 160 / $thisNodes))
        # for ((r=0; r<repeats; r++)); do
        # sbatch \
        #     --nodes=$thisNodes \
        #     --exclusive \
        #     --mem=${thisMEM}GB \
        #     ../scripts/benchmark_shmem_obs_call.sh \
        #         8 $thisSockets $thisNTPS $thisThreads $thisM $k $dem $det $flips
        # done

        # thisNodes=${shmem_ntasks16_nodes[$i]}
        # thisSockets=${shmem_ntasks16_sockets[$i]}
        # thisSPN=$((thisSockets / thisNodes))
        # thisNTPS=${shmem_ntasks16_ntps[$i]}
        # thisNTPN=$((thisSPN * thisNTPS)) # ntasks per node
        # thisMEM=$((1280))
        # for ((r=0; r<repeats; r++)); do
        # sbatch \
        #     --nodes=$thisNodes \
        #     --exclusive \
        #     --mem=${thisMEM}GB \
        #     ../scripts/benchmark_shmem_obs_call.sh \
        #         16 $thisSockets $thisNTPS $thisThreads $thisM $k $dem $det $flips
        # done
    done
    for ((i=0; i<${#shmem_quarter_threads[@]}; i++ )); do
        # Quarter subscribe sockets
        thisThreads=${shmem_quarter_threads[$i]}
        thisNTasks=${shmem_quarter_ntasks[$i]}
        thisNodes=${shmem_quarter_nodes[$i]}
        thisSockets=${shmem_quarter_sockets[$i]}
        thisNTPS=${shmem_quarter_ntps[$i]}
        thisMEM=$(($thisNTasks * 160 / $thisNodes))
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=$thisNodes \
            --exclusive \
            --mem=${thisMEM}GB \
            ../scripts/benchmark_shmem_obs_call.sh \
                $thisNTasks $thisSockets $thisNTPS $thisThreads $thisM $k $dem $det $flips
        done

        # Half subsribe sockets
        thisThreads=${shmem_half_threads[$i]}
        thisNTasks=${shmem_half_ntasks[$i]}
        thisNodes=${shmem_half_nodes[$i]}
        thisSockets=${shmem_half_sockets[$i]}
        thisNTPS=${shmem_half_ntps[$i]}
        thisMEM=$(($thisNTasks * 160 / $thisNodes))
        for ((r=0; r<repeats; r++)); do
        sbatch \
            --nodes=$thisNodes \
            --exclusive \
            --mem=${thisMEM}GB \
            ../scripts/benchmark_shmem_obs_call.sh \
                $thisNTasks $thisSockets $thisNTPS $thisThreads $thisM $k $dem $det $flips
        done
    done
done

start_serial=$(date +%s)
$serial_build predict \
    --dem $dem \
    --in $det\
    --in_format b8 \
    --out preds/preds_0.01 \
    --out_format 01 \
    --num_repeats 10 \
    >> log_0.out
end_serial=$(date +%s)
serial_time=$((end_serial - start_serial))
echo "0: $serial_time seconds" >> bench.out