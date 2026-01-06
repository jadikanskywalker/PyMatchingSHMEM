#!/usr/bin/env bash
set -euo pipefail

if ! command -v vtune >/dev/null 2>&1; then
    echo "Error: vtune command not found in PATH." >&2
    exit 1
fi

if [[ $# -lt 2 ]]; then
    cat <<'USAGE' >&2
Usage: vtune_profile.sh <result-name> <command> [args...]

Collects an Intel VTune profile for the provided command and exports
summary, top-down, and hotspots reports for later inspection.

Examples:
  scripts/vtune_profile.sh hotspot_run ./build_threads/pymatching predict --use_threads ...

Environment variables:
  VTUNE_ANALYSIS   Analysis type to collect (default: hotspots).
  VTUNE_RESULT_DIR Base directory for VTune outputs (default: run/vtune).
USAGE
    exit 1
fi

analysis=${VTUNE_ANALYSIS:-hotspots}
result_base=${VTUNE_RESULT_DIR:-run/vtune}
run_tag=$1
shift

result_dir="${result_base}/${run_tag}"
mkdir -p "${result_dir}"

vtune -collect "${analysis}" \
      -knob enable-stack-collection=true \
      -knob enable-memory-bandwidth=false \
      -r "${result_dir}/${run_tag}" \
      --app-working-dir "${PWD}" \
      -- "$@"

summary_report="${result_dir}/${run_tag}_summary.txt"
topdown_report="${result_dir}/${run_tag}_topdown.csv"
hotspots_report="${result_dir}/${run_tag}_hotspots.csv"

vtune -report summary -r "${result_dir}/${run_tag}/${run_tag}.vtune" \
      -report-output "${summary_report}"

vtune -report top-down -r "${result_dir}/${run_tag}/${run_tag}.vtune" \
      --format csv -csv-delimiter comma \
      -report-output "${topdown_report}"

vtune -report hotspots -r "${result_dir}/${run_tag}/${run_tag}.vtune" \
      --format csv -csv-delimiter comma \
      -report-output "${hotspots_report}"

cat <<EOF
VTune collection complete.
  Result dir: ${result_dir}/${run_tag}
  Summary:    ${summary_report}
  Top-down:   ${topdown_report}
  Hotspots:   ${hotspots_report}
EOF
