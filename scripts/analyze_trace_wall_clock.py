#!/usr/bin/env python3
"""Reconstruct true wall-clock phase durations from an OTF2 trace.

Computes, per shot, the wall-clock span of Local Decoding, Cross Rank Phase,
and Solution Extraction by taking min(start)/max(end) across OMP-thread
locations that touched the phase during that shot -- with zero inserted
synchronization. Shots are correlated across threads by ordinal count of
"Shot Iteration" region instances (NUM_BUFFERS_PER_UNIT=1 guarantees every
thread on a given PE works the same global shot concurrently, so no shot-id
parameter is needed in the instrumentation itself).

By default, spans are computed PER PE (merging only the threads within one
PE, never across PEs). This is always safe: a single PE's OMP threads share
one physical clock regardless of topology.

Merging locations *across* PEs into one combined span additionally requires
that the PEs share a real clock, which only holds when they run on the same
physical node (same-socket / different-socket configs). It does NOT hold
across different nodes: separate machines have independent clock epochs and
(as observed on a real cross-node run) can even progress at measurably
different effective rates, with no `--mpp=shmem` adapter here to give
Score-P any basis for cross-process clock synchronization. Merging unsynced
per-node clocks silently produces nonsensical multi-second "phase durations"
that don't correspond to anything real (caught by cross-checking against
total job wall time). Pass --merge-pes only when you know the PEs involved
share a clock domain.

Usage:
    conda run -n pymatching python scripts/analyze_trace_wall_clock.py \
        [--merge-pes] <scorep_results_dir>/traces.otf2 [<scorep_results_dir>/traces.otf2 ...]

Each anchor-file argument is one PE's trace (from scorep_wrapper.sh's
per-rank SCOREP_EXPERIMENT_DIRECTORY=scorep_results_pe<rank> layout).
"""
import argparse
import re
import statistics
import sys
from collections import defaultdict

import otf2.reader
import otf2.events

PHASE_REGIONS = ("Local Decoding", "Cross Rank Phase", "Solution Extraction")
SHOT_REGION = "Shot Iteration"
SHOT_DECODE_REGION = "Shot Decode"

# Sub-regions inside Cross Rank Phase that already isolate cross-PE communication cost as a
# pure ENTER-to-LEAVE duration on ONE PE's own clock -- no cross-PE timestamp comparison
# needed, so these remain meaningful even for cross-node runs where merging PEs is unsafe.
# "wait_until_done" (decoding_task.h) blocks on the task-status rendezvous flag; "Receiver
# Wait Until" (get_solution_from_remote_pe) blocks on the data-payload-arrived signal;
# "Sender Ctx Quiet" (send_solution_to_remote_pe) is the sender's shmem_ctx_quiet flushing
# its RMA puts before signaling done.
WAIT_REGIONS = ("wait_until_done", "Receiver Wait Until", "Sender Ctx Quiet")


def _pe_tag(anchor_path):
    """Each PE writes its own separate trace file (scorep_wrapper.sh sets a per-rank
    SCOREP_EXPERIMENT_DIRECTORY since --mpp=none treats every PE as an independent serial
    process). Location names (e.g. "Master thread") are therefore only unique *within* one
    PE's trace -- tag them with the PE's result directory name so merging multiple anchor
    paths doesn't silently collide two different PEs' "Master thread" into one bucket."""
    match = re.search(r"(scorep_results_pe\d+)", anchor_path)
    return match.group(1) if match else anchor_path


def read_intervals(anchor_path):
    """Returns {(pe_tag, location_name): [(region_name, start, end), ...]} in time order."""
    pe_tag = _pe_tag(anchor_path)
    intervals_by_location = defaultdict(list)
    with otf2.reader.open(anchor_path) as trace:
        for location, event in trace.events:
            if not isinstance(event, (otf2.events.Enter, otf2.events.Leave)):
                continue  # PROGRAM_BEGIN, THREAD_BEGIN, etc. -- not region enter/exit
            region_name = event.region.name
            key = (pe_tag, location.name)
            if isinstance(event, otf2.events.Enter):
                intervals_by_location[key].append([region_name, event.time, None])
            elif isinstance(event, otf2.events.Leave):
                # Innermost-open-region-of-this-name is the matching Enter (proper nesting
                # is guaranteed by Score-P, but different regions can interleave, so we
                # scan backward for the most recent still-open entry of the same name).
                stack = intervals_by_location[key]
                for entry in reversed(stack):
                    if entry[0] == region_name and entry[2] is None:
                        entry[2] = event.time
                        break
    return {
        loc: [(name, start, end) for name, start, end in entries]
        for loc, entries in intervals_by_location.items()
    }


def shot_windows(intervals):
    """Returns list of (shot_index, shot_start, shot_end) for one location, ordinal-ordered."""
    shots = [(start, end) for name, start, end in intervals if name == SHOT_REGION]
    shots.sort(key=lambda se: se[0])
    return [(i, start, end) for i, (start, end) in enumerate(shots)]


def union_duration(intervals):
    """Total duration covered by a list of (start, end) intervals, merging overlaps.

    Phase windows are each already a min/max across threads, so a straggler thread still
    running Local Decoding can overlap in wall-clock time with a faster thread's Cross Rank
    Phase for a different root -- summing the three phase durations directly would double
    count that overlap, so the union must be computed explicitly."""
    if not intervals:
        return 0
    ordered = sorted(intervals)
    total = 0
    cur_start, cur_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
    total += cur_end - cur_start
    return total


def phase_spans_within(intervals, region_name, window_start, window_end):
    """Returns (span_start, span_end, sum_duration, count) for instances of region_name
    whose start falls within [window_start, window_end)."""
    matches = [
        (start, end) for name, start, end in intervals
        if name == region_name and window_start <= start < window_end and end is not None
    ]
    if not matches:
        return None
    span_start = min(s for s, _ in matches)
    span_end = max(e for _, e in matches)
    sum_duration = sum(e - s for s, e in matches)
    return span_start, span_end, sum_duration, len(matches)


def _collect(anchor_paths):
    """Returns per_location, per_shot_key_phase, shot_decode_windows, per_shot_key_wait where
    per_shot_key_phase, shot_decode_windows and per_shot_key_wait are keyed by (shot_idx,
    pe_tag) -- i.e. grouping never crosses a PE boundary at this stage. Combining across PEs
    (when explicitly requested) happens only afterwards, in the caller.

    per_shot_key_wait holds, for each WAIT_REGIONS name, the *sum* of durations across all
    instances and all locations on that PE during that shot (total thread-time spent in that
    region) -- a pure single-clock quantity, safe to report even across nodes."""
    per_location = {}
    for path in anchor_paths:
        per_location.update(read_intervals(path))

    per_shot_key_phase = defaultdict(lambda: defaultdict(list))
    shot_decode_windows = defaultdict(dict)
    per_shot_key_wait = defaultdict(lambda: defaultdict(float))

    for (pe_tag, loc_name), intervals in per_location.items():
        for shot_idx, shot_start, shot_end in shot_windows(intervals):
            if shot_end is None:
                continue  # incomplete trailing shot (e.g. process exited mid-shot)
            key = (shot_idx, pe_tag)

            sd = phase_spans_within(intervals, SHOT_DECODE_REGION, shot_start, shot_end)
            if sd is not None:
                shot_decode_windows[key][loc_name] = (sd[0], sd[1])

            for phase in PHASE_REGIONS:
                result = phase_spans_within(intervals, phase, shot_start, shot_end)
                if result is not None:
                    span_start, span_end, _, _ = result
                    per_shot_key_phase[key][phase].append((loc_name, span_start, span_end))

            for wait_region in WAIT_REGIONS:
                result = phase_spans_within(intervals, wait_region, shot_start, shot_end)
                if result is not None:
                    _, _, sum_duration, _ = result
                    per_shot_key_wait[key][wait_region] += sum_duration

    return per_location, per_shot_key_phase, shot_decode_windows, per_shot_key_wait


def _phase_bounds_for_key(per_shot_key_phase, key):
    """Returns {phase: (start, end)} for one (shot_idx, pe_tag) key."""
    bounds = {}
    for phase in PHASE_REGIONS:
        entries = per_shot_key_phase[key].get(phase)
        if not entries:
            continue
        bounds[phase] = (min(s for _, s, _ in entries), max(e for _, _, e in entries))
    return bounds


def analyze_per_pe(anchor_paths, verbose=True):
    """Default, always-safe mode: spans computed within each PE separately, never merged
    across PEs. Returns (results, wait_results) where results is
    {pe_tag: {shot_idx: {phase: (start, end)}}} and wait_results is
    {pe_tag: {shot_idx: {wait_region: total_duration}}}."""
    _, per_shot_key_phase, shot_decode_windows, per_shot_key_wait = _collect(anchor_paths)

    pe_tags = sorted({pe for (_, pe) in per_shot_key_phase})
    results = {pe: {} for pe in pe_tags}
    wait_results = {pe: {} for pe in pe_tags}

    if verbose:
        print(f"{'pe':<20} {'shot':>5}  {'phase':<20} {'wall_clock_start':>18} {'wall_clock_end':>18} {'duration':>12}")
    for shot_idx, pe_tag in sorted(per_shot_key_phase, key=lambda k: (k[1], k[0])):
        key = (shot_idx, pe_tag)
        phase_bounds = _phase_bounds_for_key(per_shot_key_phase, key)
        if not phase_bounds:
            continue
        results[pe_tag][shot_idx] = phase_bounds
        wait_results[pe_tag][shot_idx] = dict(per_shot_key_wait.get(key, {}))
        if verbose:
            for phase, (start, end) in phase_bounds.items():
                print(f"{pe_tag:<20} {shot_idx:>5}  {phase:<20} {start:>18} {end:>18} {end - start:>12}")
            for wait_region, total in wait_results[pe_tag][shot_idx].items():
                print(f"{pe_tag:<20} {shot_idx:>5}  {'  [wait] ' + wait_region:<20} {'':>18} {'':>18} {total:>12.0f}")

            sd_windows = shot_decode_windows.get(key, {})
            if sd_windows:
                sd_start = min(s for s, _ in sd_windows.values())
                sd_end = max(e for _, e in sd_windows.values())
                covered = union_duration(list(phase_bounds.values()))
                residual = (sd_end - sd_start) - covered
                print(f"{pe_tag:<20} {shot_idx:>5}  {'(residual/gap)':<20} {'':>18} {'':>18} {residual:>12}")

    return results, wait_results


def analyze_merged(anchor_paths):
    """Opt-in mode: merges spans across ALL given PEs into one wall-clock window per shot.
    Only valid when every PE involved shares a real clock (i.e. same physical node) --
    see module docstring. Returns per_shot_phase, shot_decode_windows (both keyed by
    shot_idx only, pooling every PE together)."""
    print("WARNING: --merge-pes assumes all PEs share a real clock (same physical node).",
          file=sys.stderr)
    print("         Do not trust this mode's output for cross-node configurations.",
          file=sys.stderr)

    _, per_shot_key_phase, shot_decode_windows_by_key, _ = _collect(anchor_paths)

    # Pool every PE together under just shot_idx.
    per_shot_phase = defaultdict(lambda: defaultdict(list))
    shot_decode_windows = defaultdict(dict)
    for (shot_idx, pe_tag), phases in per_shot_key_phase.items():
        for phase, entries in phases.items():
            per_shot_phase[shot_idx][phase].extend(entries)
    for (shot_idx, pe_tag), locs in shot_decode_windows_by_key.items():
        shot_decode_windows[shot_idx].update({f"{pe_tag}:{loc}": v for loc, v in locs.items()})

    print(f"{'shot':>5}  {'phase':<20} {'wall_clock_start':>18} {'wall_clock_end':>18} {'duration':>12}")
    for shot_idx in sorted(per_shot_phase):
        phase_bounds = {}
        for phase in PHASE_REGIONS:
            entries = per_shot_phase[shot_idx].get(phase)
            if not entries:
                continue
            start = min(s for _, s, _ in entries)
            end = max(e for _, _, e in entries)
            phase_bounds[phase] = (start, end)
            print(f"{shot_idx:>5}  {phase:<20} {start:>18} {end:>18} {end - start:>12}")

        sd_windows = shot_decode_windows.get(shot_idx, {})
        if sd_windows and phase_bounds:
            sd_start = min(s for s, _ in sd_windows.values())
            sd_end = max(e for _, e in sd_windows.values())
            covered = union_duration(list(phase_bounds.values()))
            residual = (sd_end - sd_start) - covered
            print(f"{shot_idx:>5}  {'(residual/gap)':<20} {'':>18} {'':>18} {residual:>12}")

    return per_shot_phase, shot_decode_windows


def summarize_per_pe(results, wait_results, warmup_shots=2):
    """Prints median/mean/max phase duration per PE, excluding the first `warmup_shots`.

    Wait-region figures are per-shot *totals* summed across all threads on that PE (thread-
    time spent waiting, not a wall-clock span) -- see WAIT_REGIONS docstring."""
    for pe_tag in sorted(results):
        print(f"=== {pe_tag} ({len(results[pe_tag])} shots total) ===")
        for phase in PHASE_REGIONS:
            durs = []
            for shot_idx, bounds in results[pe_tag].items():
                if shot_idx < warmup_shots:
                    continue
                if phase in bounds:
                    start, end = bounds[phase]
                    durs.append((end - start) / 1e6)  # ms
            if durs:
                print(f"  {phase:<20} n={len(durs):4d}  median={statistics.median(durs):9.3f}ms  "
                      f"mean={statistics.mean(durs):9.3f}ms  max={max(durs):9.3f}ms")
        for wait_region in WAIT_REGIONS:
            durs = []
            for shot_idx, waits in wait_results.get(pe_tag, {}).items():
                if shot_idx < warmup_shots:
                    continue
                if wait_region in waits:
                    durs.append(waits[wait_region] / 1e6)  # ms
            if durs:
                print(f"  [wait] {wait_region:<20} n={len(durs):4d}  median={statistics.median(durs):9.3f}ms  "
                      f"mean={statistics.mean(durs):9.3f}ms  max={max(durs):9.3f}ms")
        print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("traces", nargs="+", help="one or more scorep_results_pe*/traces.otf2 paths")
    parser.add_argument("--merge-pes", action="store_true",
                         help="merge spans across all given PEs (only valid when they share a clock domain, e.g. same node)")
    parser.add_argument("--summary-only", action="store_true",
                         help="skip the per-shot table, print only the median/mean/max summary")
    args = parser.parse_args()

    if args.merge_pes:
        analyze_merged(args.traces)
    else:
        results, wait_results = analyze_per_pe(args.traces, verbose=not args.summary_only)
        summarize_per_pe(results, wait_results)
