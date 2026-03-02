import MDAnalysis as mda
import numpy as np
import sys
import os
from pathlib import Path

def main():
    # 1. Check if the user provided a input_path
    if len(sys.argv) < 2:
        print("Usage: python3 converter.py <input_file> (optional)<coords_output> (optional)<meta_output>")
        sys.exit(1)

    # Generate output names based on input name
    input_path=sys.argv[1]
    bin_output = sys.argv[2] if len(sys.argv)>2 else f"temp_coords.bin"
    meta_output = sys.argv[3] if len(sys.argv)>3 else f"temp_metadata.txt"
    try:
        # 2. Load the Universe
        u = mda.Universe(input_path)
        positions = u.atoms.positions.astype(np.float64)
#        print(positions[:20,:]) # for checking

        # 3. Save binary coordinates
        positions.tofile(bin_output)

        # 4. Save metadata (N atoms)
        ext_no_dot = Path(input_path).suffix[1:]
        with open(meta_output, "w") as f:
            f.write(f"{str(len(u.atoms))}\n")
            f.write(f"{ext_no_dot}\n")
            if ext_no_dot=="gro":
              dimensions = u.dimensions[:3] # [Lx, Ly, Lz] in Angstroms
              f.write(f"{len(u.atoms)} {dimensions[0]} {dimensions[1]} {dimensions[2]}")

        print(f"Successfully converted {len(u.atoms)} atoms.")
        print(f"Output: {bin_output}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()