#!/bin/bash
#SBATCH --job-name=profile_threads
#SBATCH --partition=zen4
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --sockets=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=64G

if [ $# -lt 2 ]
  then
    echo "Args: [nthreads] [M]"
    exit 1
else
    nthreads=$1
    M=$2
fi

source ~/.bash_profile
conda activate pymatching

FG_DIR=~/sw/el9-x86_64/FlameGraph
threads_build=~/PyMatchingSHMEM/build_threads_release/pymatching

export OMP_NUM_THREADS=$nthreads
export OMP_PLACES=cores
export OMP_PROC_BIND=close

tag=threads${nthreads}

~/PyMatchingSHMEM/scripts/vtune_profile.sh $tag \
    $threads_build predict \
    --dem ../error_model.dem \
    --in ../detection_events.b8 \
    --in_format b8 \
    --out predicted_obs_flips.01 \
    --out_format 01 \
    --rounds_per_partition $M \
    --use_threads
if [ $? -ne 0 ]; then
    echo "vtune_profile.sh failed for $tag" >&2
    exit 1
fi

echo "correct predictions:"
paste -d " " predicted_obs_flips.01 ../actual_obs_flips.01 | grep "1 1\|0 0" | wc -l
echo "wrong predictions:"
paste -d " " predicted_obs_flips.01 ../actual_obs_flips.01 | grep "0 1\|1 0" | wc -l

topdown_csv=run/vtune/${tag}/${tag}_topdown.csv
sed -i '1{/^war:/d}' "$topdown_csv"
perl $FG_DIR/stackcollapse-vtune.pl "$topdown_csv" \
    | perl $FG_DIR/flamegraph.pl --width 1400 --fontsize 10 \
    > run/vtune/${tag}/${tag}_flamegraph.svg

echo "Flamegraph: $(pwd)/run/vtune/${tag}/${tag}_flamegraph.svg"

# Source-line breakdown of decode_shots's self time -- with -g in the release
# build, this splits out inlined helpers (wait_until_decode_done etc.) that
# the function-level hotspots report bundles into one flat bucket.
sourceline_report=run/vtune/${tag}/${tag}_decode_shots_sourceline.csv
vtune -report hotspots -r run/vtune/${tag}/${tag}/${tag}.vtune \
    -group-by "source-line" \
    -filter "function=pm::DecodingUnit::decode_shots._omp_fn.0" \
    --format csv -csv-delimiter comma \
    -report-output "$sourceline_report"
if [ $? -ne 0 ]; then
    echo "source-line report failed for $tag" >&2
fi
