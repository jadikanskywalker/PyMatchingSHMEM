#!/usr/bin/env python3
"""Scan out_parallel/p*/t*.out debug traces for a task solved more than once in the
same shot -- either by the same thread twice or by two different threads.
Also reports the shot context around any 'ERROR:' lines found.
"""
import re
import sys
import glob
from collections import defaultdict

SHOT_RE = re.compile(r"Starting shot (\d+), buffer round (\d+), buffer_id (\d+)")
SOLVING_RE = re.compile(r"Thread (\d+) solving ([pf])(\d+)")
ERROR_RE = re.compile(r"^ERROR:")


def scan_dir(out_parallel_dir):
    # (pe, shot_id, kind, part) -> list of (thread_file, thread_id, line_no)
    solves = defaultdict(list)
    errors = []  # (pe, thread_file, shot_id, line_no, text)

    for pe_dir in sorted(glob.glob(f"{out_parallel_dir}/p*")):
        pe = pe_dir.split("/")[-1]
        for t_file in sorted(glob.glob(f"{pe_dir}/t*.out")):
            shot_id = None
            with open(t_file, errors="replace") as f:
                for line_no, line in enumerate(f, 1):
                    m = SHOT_RE.search(line)
                    if m:
                        shot_id = int(m.group(1))
                        continue
                    m = SOLVING_RE.search(line)
                    if m:
                        tid, kind, part = int(m.group(1)), m.group(2), int(m.group(3))
                        solves[(pe, shot_id, kind, part)].append((t_file, tid, line_no))
                        continue
                    if ERROR_RE.search(line):
                        errors.append((pe, t_file, shot_id, line_no, line.rstrip()))

    return solves, errors


def main():
    if len(sys.argv) != 2:
        print("Usage: check_double_solve.py <out_parallel_dir>")
        sys.exit(1)
    out_parallel_dir = sys.argv[1]
    solves, errors = scan_dir(out_parallel_dir)

    dupes = {k: v for k, v in solves.items() if len(v) > 1}
    print(f"Total distinct (pe, shot, part) solve-events: {len(solves)}")
    print(f"DUPLICATE solves found: {len(dupes)}")
    for (pe, shot_id, kind, part), occurrences in sorted(dupes.items(), key=lambda kv: (kv[0][1] or -1, kv[0][3])):
        print(f"  {pe} shot={shot_id} {kind}{part} solved {len(occurrences)}x:")
        for t_file, tid, line_no in occurrences:
            print(f"    {t_file}:{line_no} (thread {tid})")

    print()
    print(f"ERROR lines found: {len(errors)}")
    # Report the shot(s) in which errors occurred, and whether any dupes overlap that shot
    error_shots = defaultdict(int)
    for pe, t_file, shot_id, line_no, text in errors:
        error_shots[(pe, shot_id)] += 1
    for (pe, shot_id), count in sorted(error_shots.items(), key=lambda kv: (kv[0][1] or -1)):
        print(f"  {pe} shot={shot_id}: {count} error line(s)")
        overlap = [k for k in dupes if k[0] == pe and k[1] == shot_id]
        if overlap:
            print(f"    -> {len(overlap)} duplicate solve(s) also in this shot: {overlap}")

    if not dupes and not errors:
        print("No duplicates or errors found.")


if __name__ == "__main__":
    main()
