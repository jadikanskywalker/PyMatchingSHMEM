import re
import os
import glob
import sys

def parse_log_file(filepath):
    """Returns max decoding time found in the file."""
    times = []
    with open(filepath, 'r') as f:
        for line in f:
            match = re.search(r'Decoding time:\s*([0-9.]+)\s*ms', line)
            if match:
                times.append(float(match.group(1)))
    if not times:
        return None
    return max(times)

def main():
    directory = 'bench_obs_001'
    if len(sys.argv) > 1:
        directory = sys.argv[1]
        
    results = []
    
    for filepath in glob.glob(os.path.join(directory, 'log_*.out')):
        filename = os.path.basename(filepath)
        decoding_time = parse_log_file(filepath)
        if decoding_time is None:
            continue
            
        m_val, ntasks_val, sockets_val, ntps_val, nthreads_val = "", 1, 1, 1, 1

        if filename == 'log_0.out':
            m_val, ntasks_val, sockets_val, ntps_val, nthreads_val = "", 1, 1, 1, 1
        else:
            match = re.match(r'log_M(\d+)_ntasks(\d+)_sockets(\d+)_ntps(\d+)_nthreads(\d+)_k\d+', filename)
            if match:
                m_val        = int(match.group(1))
                ntasks_val   = int(match.group(2))
                sockets_val  = int(match.group(3))
                ntps_val     = int(match.group(4))
                nthreads_val = int(match.group(5))
            else:
                continue

        nodes_val = max(sockets_val / 2, 1)

        results.append({
            'M':             m_val,
            'ntasks':        ntasks_val,
            'sockets':       sockets_val,
            'nodes':         nodes_val,
            'ntps':          ntps_val,
            'nthreads':      nthreads_val,
            'Decoding Time': decoding_time
        })

    def sort_key(row):
        m = 0 if row['M'] == "" else row['M']
        return (m, row['ntasks'], row['sockets'], row['nthreads'])

    results.sort(key=sort_key)

    output_file = os.path.join(directory, 'results.csv')
    with open(output_file, 'w') as f:
        f.write("M,ntasks,sockets,nodes,ntsp,nthreads,Decoding Time\n")

        prev_m = None
        for r in results:
            display_m = r['M'] if r['M'] != prev_m else ""
            prev_m = r['M']
            f.write(f"{display_m},{r['ntasks']},{r['sockets']},{r['nodes']},{r['ntps']},{r['nthreads']},{r['Decoding Time']}\n")
            
    print(f"Results written to {output_file}")

if __name__ == '__main__':
    main()
