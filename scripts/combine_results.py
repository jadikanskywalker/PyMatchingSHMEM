import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(description="Combine PyMatching parallel output files by XORing them.")
    parser.add_argument("base_filename", help="Base filename of the output (e.g. out.01)")
    parser.add_argument("num_pes", type=int, help="Number of PEs/files to combine")
    parser.add_argument("--format", choices=["01", "b8"], default="01", help="Format of the output files (01 or b8)")
    
    args = parser.parse_args()
    
    base_name = args.base_filename
    # Handle extension if present in base_name logic of C++ code
    # The C++ code inserts _pe<pid> before the last dot.
    
    parts = os.path.splitext(base_name)
    prefix = parts[0]
    ext = parts[1]
    
    filenames = []
    for pid in range(args.num_pes):
        fn = f"{prefix}_pe{pid}{ext}"
        if not os.path.exists(fn):
            print(f"Error: File {fn} not found.", file=sys.stderr)
            sys.exit(1)
        filenames.append(fn)

    print(f"Combining {len(filenames)} files: {filenames} -> {base_name}")

    if args.format == "01":
        combine_01(filenames, base_name)
    elif args.format == "b8":
        combine_b8(filenames, base_name)

def combine_01(filenames, out_filename):
    files = [open(fn, "r") for fn in filenames]
    with open(out_filename, "w") as out_f:
        while True:
            lines = [f.readline() for f in files]
            if not lines[0]:
                break # EOF
            
            # Check all lines have same length
            line_len = len(lines[0])
            for i, line in enumerate(lines):
                 if len(line) != line_len:
                      print(f"Warning: Line length mismatch in file {filenames[i]}", file=sys.stderr)

            # XOR
            # Assuming lines end with newline
            result_chars = []
            for i in range(len(lines[0].strip())):
                val = 0
                for line in lines:
                    if i < len(line) and line[i] in '01':
                        val ^= int(line[i])
                result_chars.append(str(val))
            
            out_f.write("".join(result_chars) + "\n")
            
    for f in files:
        f.close()

def combine_b8(filenames, out_filename):
    # b8 format: bits packed into bytes.
    # We can read chunks and XOR bytes.
    chunk_size = 4096
    files = [open(fn, "rb") for fn in filenames]
    with open(out_filename, "wb") as out_f:
        while True:
            chunks = [f.read(chunk_size) for f in files]
            if not chunks[0]:
                break
            
            # Use bytearray for mutability or just create new bytes
            # XOR byte by byte
            # Assuming all read same amount
            res = bytearray(len(chunks[0]))
            for i in range(len(chunks[0])):
                val = 0
                for chunk in chunks:
                    if i < len(chunk):
                        val ^= chunk[i]
                res[i] = val
            out_f.write(res)
            
    for f in files:
        f.close()

if __name__ == "__main__":
    main()
