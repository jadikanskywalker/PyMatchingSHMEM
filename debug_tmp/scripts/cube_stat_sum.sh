#!/bin/bash
# For a bench_36obs_crt_scaling-style dir: runs cube_stat per PE (cube_stat_pe<N>.out, same as
# manual runs this session), then sums NumberOfCalls/ExclusiveTime/InclusiveTime per region name
# across all PEs into cube_stat_summary.out (sorted by total ExclusiveTime, descending).
set -e
if [ $# -ne 1 ]; then
    echo "Args: [dir containing scorep_results_pe*/profile.cubex]"
    exit 1
fi
dir=$1

for cubex in "$dir"/scorep_results_pe*/profile.cubex; do
    pe=$(basename "$(dirname "$cubex")" | sed 's/scorep_results_//')
    cube_stat -t inf -p "$cubex" > "$dir/cube_stat_${pe}.out"
done

{
    printf "%-28s %13s %13s %13s\n" "cube::Region" "NumberOfCalls" "ExclusiveTime" "InclusiveTime"
    awk '
    FNR==1 { next }  # skip header line in each per-PE file
    {
        n = NF
        incl = $n; excl = $(n-1); calls = $(n-2)
        name = $1
        for (i = 2; i <= n-3; i++) name = name " " $i
        total_calls[name] += calls
        total_excl[name]  += excl
        total_incl[name]  += incl
    }
    END {
        for (name in total_excl)
            print total_excl[name], total_calls[name], total_incl[name], name
    }' "$dir"/cube_stat_pe*.out \
    | sort -t' ' -k1,1 -rn \
    | awk '{ name=$4; for (i=5; i<=NF; i++) name = name " " $i
             printf "%-28s %13.0f %13.6f %13.6f\n", name, $2, $1, $3 }'
} > "$dir/cube_stat_summary.out"

echo "Wrote $dir/cube_stat_summary.out"
