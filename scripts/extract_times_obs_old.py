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
            
        m_val, nodes_val, ppn_val, nthreads_val = "", 1, 1, 1
        
        if filename == 'log_0.out':
            m_val = ""
            nodes_val = 1
            pps_val = 1
            np_val = 1
            nthreads_val = 1
        else:
            # Check for shmem format: log_M32_n4_ppn2_32threads_k1.out
            # Check for thread format: log_M32_32threads.out
            
            shmem_match = re.match(r'log_M(\d+)_n(\d+)_pps(\d+)_(\d+)_k(\d+)', filename)
            # thread_match = re.match(r'log_M(\d+)_(\d+)threads\.out', filename)
            
            if shmem_match:
                m_val = int(shmem_match.group(1))
                np_val = int(shmem_match.group(2))
                pps_val = int(shmem_match.group(3))
                nthreads_val = int(shmem_match.group(4))
                nodes_val = max(np_val // (pps_val*2), 1)
            # elif thread_match:
            #     m_val = int(thread_match.group(1))
            #     nthreads_val = int(thread_match.group(2))
            #     nodes_val = 1
            #     ppn_val = 1
            else:
                continue
                
        results.append({
            'M': m_val,
            'cpus_total': nthreads_val*np_val,
            'nodes': nodes_val,
            'ntasks': np_val,
            'cpus_per_task': nthreads_val,
            'Decoding Time': decoding_time
        })
        
    # Sort the results according to the table logic
    def sort_key(row):
        m = 0 if row['M'] == "" else row['M']
        return (m, row['cpus_total'], row['ntasks'], row['nodes'], row['cpus_per_task'])

    results.sort(key=sort_key)
    
    output_file = os.path.join(directory, 'results.csv')
    with open(output_file, 'w') as f:
        f.write("M,nodes,ntasks,cpus_per_task,cpus_total,Decoding Time\n")
        
        prev_m = None
        for r in results:
            display_m = r['M'] if r['M'] != prev_m else ""
            prev_m = r['M']
            f.write(f"{display_m},{r['nodes']},{r['ntasks']},{r['cpus_per_task']},{r['cpus_total']},{r['Decoding Time']}\n")
            
    print(f"Results written to {output_file}")

if __name__ == '__main__':
    main()
