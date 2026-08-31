#!/bin/bash
#SBATCH --job-name=regen_flamegraphs
#SBATCH --output=debug_tmp/out/regen_flamegraphs-%j.out
#SBATCH --error=debug_tmp/out/regen_flamegraphs-%j.err
#SBATCH --partition=zen4
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

source ~/.bash_profile
conda activate pymatching

FG_DIR=~/sw/el9-x86_64/FlameGraph
base=~/PyMatchingSHMEM/run_d21_profile

for t in 16 32 64 128 256; do
    dir="${base}/log_${t}threads/run/vtune/threads${t}"
    tag="threads${t}"
    vtune -report top-down -r "${dir}/${tag}/${tag}.vtune" \
          -call-stack-mode all -column="CPU Time:Self","Module" -filter "Function Stack" \
          --format csv -csv-delimiter comma \
          -report-output "${dir}/${tag}_topdown.csv"
    if [ $? -ne 0 ]; then
        echo "vtune report failed for ${tag}" >&2
        continue
    fi

    sed -i '1{/^war:/d}' "${dir}/${tag}_topdown.csv"
    perl $FG_DIR/stackcollapse-vtune.pl "${dir}/${tag}_topdown.csv" \
        | perl $FG_DIR/flamegraph.pl --width 1400 --fontsize 10 \
        > "${dir}/${tag}_flamegraph.svg"
    echo "Regenerated ${dir}/${tag}_flamegraph.svg"
done
