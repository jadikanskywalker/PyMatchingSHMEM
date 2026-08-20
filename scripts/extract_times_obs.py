import glob
import os
import re
import sys

# Parses the per-PE "Average decoding time: X.XXms" summary line each SHMEM run prints into its own
# log_shmem_pe<pid>.out (one config subdir per run, out_M<M>_ntasks<...>_sockets<...>_ntps<...>_
# nthreads<...>_k<...>_<jobid>/). Takes the MAX across a config's own PEs (the slowest PE is the
# real wall-clock bottleneck). A config with zero matching PE logs (construction-time "not enough
# threads" throw, or a launch/ppr placement failure) is skipped entirely -- there's no timing data
# to report, not a zero. Older versions of this script (and its sibling extract_times.py) scanned
# for a "Decoding time: X ms" print that doesn't exist in the current SHMEM pipeline; wall-clock
# `date +%s` numbers in bench.out looked like a substitute but silently include failed runs (e.g. a
# construction-throw config still burns real wall-clock time before exiting) -- this goes back to
# the actual per-PE decode-time print instead, which only exists when a run genuinely decoded.

CONFIG_DIR_RE = re.compile(
    r'^out_M(\d+)_ntasks(\d+)_sockets(\d+)_ntps(\d+)_nthreads(\d+)_k\d+(?:_\d+)?$')
AVG_RE = re.compile(r'Average decoding time:\s*([0-9.]+)\s*ms')
SERIAL_RE = re.compile(r'Decoding time:\s*([0-9.]+)\s*ms')


def max_avg_decoding_time(config_dir):
    times = []
    for pe_log in glob.glob(os.path.join(config_dir, 'log_shmem_pe*.out')):
        with open(pe_log, 'r') as f:
            for line in f:
                m = AVG_RE.search(line)
                if m:
                    times.append(float(m.group(1)))
    return max(times) if times else None


def serial_decoding_time(filepath):
    if not os.path.exists(filepath):
        return None
    time = None
    with open(filepath, 'r') as f:
        for line in f:
            m = SERIAL_RE.search(line)
            if m:
                time = float(m.group(1))
    return time


def main():
    directory = 'bench_obs_001'
    if len(sys.argv) > 1:
        directory = sys.argv[1]

    results = []

    serial_time = serial_decoding_time(os.path.join(directory, 'log_0.out'))
    if serial_time is not None:
        results.append({
            'M': "", 'ntasks': 1, 'sockets': 1, 'nodes': 1, 'ntps': 1, 'nthreads': 1,
            'Decoding Time (ms)': serial_time,
        })

    for entry in sorted(os.listdir(directory)):
        m = CONFIG_DIR_RE.match(entry)
        if not m:
            continue
        m_val, ntasks, sockets, ntps, nthreads = (int(x) for x in m.groups())
        config_dir = os.path.join(directory, entry)
        decoding_time = max_avg_decoding_time(config_dir)
        if decoding_time is None:
            continue
        results.append({
            'M': m_val,
            'ntasks': ntasks,
            'sockets': sockets,
            'nodes': max(sockets / 2, 1),
            'ntps': ntps,
            'nthreads': nthreads,
            'Decoding Time (ms)': decoding_time,
        })

    def sort_key(row):
        m = 0 if row['M'] == "" else row['M']
        return (m, row['ntasks'], row['sockets'], row['nthreads'])

    results.sort(key=sort_key)

    output_file = os.path.join(directory, 'results.csv')
    with open(output_file, 'w') as f:
        f.write("M,ntasks,sockets,nodes,ntps,nthreads,Decoding Time (ms)\n")
        prev_m = None
        for r in results:
            display_m = r['M'] if r['M'] != prev_m else ""
            prev_m = r['M']
            f.write(f"{display_m},{r['ntasks']},{r['sockets']},{r['nodes']},{r['ntps']},"
                    f"{r['nthreads']},{r['Decoding Time (ms)']}\n")

    print(f"Results written to {output_file} ({len(results)} rows)")


if __name__ == '__main__':
    main()
