#!/usr/bin/env python3
"""
Aggregate Score-P cube_stat CSV profiles across all PEs for a given run name.

Usage:
    python3 scripts/aggregate_profiles.py <name>

Example:
    python3 scripts/aggregate_profiles.py 4obs

For each of run_same_socket_<name>, run_diff_socket_<name>, run_diff_node_<name>,
finds all scorep_results_pe*/profile.cubex, runs cube_stat -t 20 on each, sums
Number of Calls / Exclusive Time / Inclusive Time per region across all PEs, and
writes aggregated_profile.csv to the run directory.
"""

import sys
import csv
import subprocess
import io
from collections import defaultdict
from pathlib import Path


CONFIGS = ["same_socket", "diff_socket", "diff_node", "single_pe"]
CUBE_STAT_THREADS = 20


def run_cube_stat(cubex_path: Path) -> list[dict]:
    result = subprocess.run(
        ["cube_stat", "-t", str(CUBE_STAT_THREADS), str(cubex_path)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"cube_stat failed on {cubex_path}:\n{result.stderr.strip()}")
    reader = csv.DictReader(io.StringIO(result.stdout))
    rows = []
    for row in reader:
        rows.append({
            "region": row["cube::Region"].strip('"'),
            "calls": int(row["Number of Calls"]),
            "excl": float(row["Exclusive Time"]),
            "incl": float(row["Inclusive Time"]),
        })
    return rows


def aggregate_pe_profiles(run_dir: Path) -> list[dict]:
    pe_dirs = sorted(run_dir.glob("scorep_results_pe*/profile.cubex"))
    if not pe_dirs:
        raise FileNotFoundError(f"No scorep_results_pe*/profile.cubex found in {run_dir}")

    totals: dict[str, dict] = defaultdict(lambda: {"calls": 0, "excl": 0.0, "incl": 0.0})
    for cubex in pe_dirs:
        print(f"  Reading {cubex.relative_to(run_dir.parent)} ...", flush=True)
        for row in run_cube_stat(cubex):
            totals[row["region"]]["calls"] += row["calls"]
            totals[row["region"]]["excl"]  += row["excl"]
            totals[row["region"]]["incl"]  += row["incl"]

    # Sort by descending inclusive time (same ordering as cube_stat default)
    return sorted(
        [{"region": r, **v} for r, v in totals.items()],
        key=lambda x: x["incl"], reverse=True
    )


def write_csv(rows: list[dict], out_path: Path) -> None:
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(["cube::Region", "Number of Calls", "Exclusive Time", "Inclusive Time"])
        for row in rows:
            writer.writerow([
                row["region"],
                row["calls"],
                round(row["excl"], 6),
                round(row["incl"], 6),
            ])
    print(f"  -> {out_path}")


def main():
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    name = sys.argv[1]
    repo_root = Path(__file__).resolve().parent.parent

    found_any = False
    for config in CONFIGS:
        run_dir = repo_root / f"run_{config}_{name}"
        print(f"\n[{config}]  {run_dir.name}")
        if not run_dir.is_dir():
            print(f"  Skipping: directory not found")
            continue
        try:
            rows = aggregate_pe_profiles(run_dir)
            out_path = run_dir / "aggregated_profile.csv"
            write_csv(rows, out_path)
            found_any = True
        except FileNotFoundError as e:
            print(f"  Skipping: {e}")
        except RuntimeError as e:
            print(f"  Error: {e}", file=sys.stderr)

    if not found_any:
        print(f"\nNo run directories found for name '{name}'.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
