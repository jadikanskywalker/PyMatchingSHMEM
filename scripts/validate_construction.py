#!/usr/bin/env python3
"""Parse the DEBUG construction dump (log_shmem_pe*.out) and validate the task graph:
parent/child consistency, root count, duplicate part ids, degenerate seam windows,
and how many distinct Tasks each special task (seam/CRT) is attached to.
"""
import re
import sys
from collections import defaultdict

PART_RE = re.compile(r"^--part: (\S+)\s+vb_left: (-?\d+)\s+vb_right: (-?\d+)\s+vb_marker: (-?\d+)\s+solver: (\S+)")
FLAGS_RE = re.compile(r"is_extraction_unit_root: (\d)\s+is_extraction_unit_connector: (\d)\s+defer_division: (\d)")
CHILDREN_RE = re.compile(r"left_child: (0x[0-9a-fA-F]+)(?:\([^)]*\))?\s+me: (0x[0-9a-fA-F]+)\s+right_child: (0x[0-9a-fA-F]+)(?:\([^)]*\))?")
PARENT_RE = re.compile(r"parent: f?(0x[0-9a-fA-F]+)")
SPECIAL_HDR_RE = re.compile(r"special_tasks \((\d+)\):(.*)")
SPECIAL_ENTRY_RE = re.compile(r"(0x[0-9a-fA-F]+)\((seam|crt) vb=(-?\d+) vb_left=(-?\d+) vb_right=(-?\d+)\)")
NUM_ROOTS_RE = re.compile(r"num_task_roots=(\d+)")


def parse(path):
    tasks = {}  # me_addr -> dict
    special_owner = defaultdict(list)  # special_task_addr -> list of owning task addrs
    special_info = {}  # special_task_addr -> (kind, vb, vb_left, vb_right)
    num_task_roots = None

    with open(path, errors="replace") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        m = NUM_ROOTS_RE.search(line)
        if m:
            num_task_roots = int(m.group(1))
        m = PART_RE.match(line)
        if m:
            part_name, vb_left, vb_right, vb_marker, solver = m.groups()
            flags_m = FLAGS_RE.search(lines[i + 1])
            children_m = CHILDREN_RE.search(lines[i + 2])
            parent_m = PARENT_RE.search(lines[i + 3])
            special_m = SPECIAL_HDR_RE.search(lines[i + 4])
            if not (flags_m and children_m and parent_m and special_m):
                print(f"WARNING: malformed entry near line {i+1} for {part_name}")
                i += 1
                continue
            me_addr = children_m.group(2)
            entry = {
                "part_name": part_name,
                "vb_left": int(vb_left),
                "vb_right": int(vb_right),
                "vb_marker": int(vb_marker),
                "is_root": flags_m.group(1) == "1",
                "is_connector": flags_m.group(2) == "1",
                "defer_division": flags_m.group(3) == "1",
                "left_child": children_m.group(1),
                "right_child": children_m.group(3),
                "parent": parent_m.group(1),
                "line": i + 1,
            }
            if me_addr in tasks:
                print(f"DUPLICATE me_addr {me_addr} at line {i+1} (first seen for {tasks[me_addr]['part_name']})")
            tasks[me_addr] = entry

            n_special = int(special_m.group(1))
            rest = special_m.group(2)
            found = SPECIAL_ENTRY_RE.findall(rest)
            if len(found) != n_special:
                print(f"WARNING: special_tasks count mismatch for {part_name} (line {i+1}): header says {n_special}, parsed {len(found)}")
            for addr, kind, vb, vbl, vbr in found:
                special_owner[addr].append(me_addr)
                special_info[addr] = (kind, int(vb), int(vbl), int(vbr))
            i += 5
            continue
        i += 1

    return tasks, special_owner, special_info, num_task_roots


def main():
    if len(sys.argv) != 2:
        print("Usage: validate_construction.py <log_shmem_peN.out>")
        sys.exit(1)
    path = sys.argv[1]
    tasks, special_owner, special_info, num_task_roots = parse(path)
    print(f"Parsed {len(tasks)} tasks. num_task_roots (reported) = {num_task_roots}")

    # 1. Parent/child pointer resolution + bidirectional consistency
    dangling = 0
    inconsistent = 0
    roots = []
    for addr, t in tasks.items():
        for side, child_addr in [("left", t["left_child"]), ("right", t["right_child"])]:
            if child_addr != "0x0" and child_addr not in tasks:
                print(f"DANGLING {side}_child {child_addr} on {t['part_name']} (line {t['line']})")
                dangling += 1
        if t["parent"] == "0x0":
            roots.append(addr)
        else:
            if t["parent"] not in tasks:
                print(f"DANGLING parent {t['parent']} on {t['part_name']} (line {t['line']})")
                dangling += 1
            else:
                p = tasks[t["parent"]]
                if p["left_child"] != addr and p["right_child"] != addr:
                    print(f"INCONSISTENT: {t['part_name']} (line {t['line']}) claims parent {t['parent']} "
                          f"({p['part_name']}), but that parent's children are {p['left_child']}/{p['right_child']}, not {addr}")
                    inconsistent += 1

    print(f"Dangling pointers: {dangling}")
    print(f"Parent/child inconsistencies: {inconsistent}")
    print(f"Actual root count (parent==0x0): {len(roots)}", "MATCHES" if len(roots) == num_task_roots else "MISMATCH vs num_task_roots!")

    # 2. Duplicate part numbers
    part_seen = defaultdict(list)
    for addr, t in tasks.items():
        part_seen[t["part_name"]].append(addr)
    dup_parts = {k: v for k, v in part_seen.items() if len(v) > 1}
    print(f"Duplicate part names: {len(dup_parts)}")
    for k, v in list(dup_parts.items())[:10]:
        print(f"  {k}: {v}")

    # 3. Degenerate seam/crt windows
    degenerate = [(addr, info) for addr, info in special_info.items() if info[3] <= info[2]]
    print(f"Degenerate (vb_right <= vb_left) special tasks: {len(degenerate)}")
    for addr, (kind, vb, vbl, vbr) in degenerate[:10]:
        print(f"  {addr} {kind} vb={vb} vb_left={vbl} vb_right={vbr} owners={special_owner[addr]}")

    # 4. Special task owner-count distribution (expect 2 per seam in the common case)
    owner_counts = defaultdict(int)
    for addr, owners in special_owner.items():
        owner_counts[len(owners)] += 1
    print(f"Special-task owner-count distribution: {dict(owner_counts)}")
    weird = {addr: owners for addr, owners in special_owner.items() if len(owners) == 1}
    if weird:
        print(f"Special tasks attached to only ONE task (never resolvable): {len(weird)}")
        for addr, owners in list(weird.items())[:10]:
            kind, vb, vbl, vbr = special_info[addr]
            print(f"  {addr} {kind} vb={vb} owner={owners}")

    # 5. Leaf partition coverage
    leaf_nums = []
    for addr, t in tasks.items():
        if t["part_name"].startswith("p"):
            leaf_nums.append(int(t["part_name"][1:]))
    leaf_nums.sort()
    expected = list(range(len(leaf_nums)))
    if leaf_nums != expected:
        missing = sorted(set(expected) - set(leaf_nums))
        extra = sorted(set(leaf_nums) - set(expected))
        print(f"LEAF COVERAGE MISMATCH: {len(leaf_nums)} leaves found, expected {len(expected)} contiguous from 0.")
        print(f"  missing: {missing[:20]}")
        print(f"  extra/duplicate: {extra[:20]}")
    else:
        print(f"Leaf coverage OK: {len(leaf_nums)} leaves, contiguous 0..{len(leaf_nums)-1}")


if __name__ == "__main__":
    main()
