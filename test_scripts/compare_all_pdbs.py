import numpy as np
import re

# ============================================================
# Load coordinates from PDB (robust)
# ============================================================
def load_pdb(path):
    coords = []
    atom_pattern = re.compile(r'(ATOM|HETATM)')
    try:
        with open(path, 'r') as f:
            for line in f:
                if not atom_pattern.search(line):
                    continue

                # Extract using fixed columns
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except:
                    # fallback: split-based extraction
                    parts = line.split()
                    nums = [p for p in parts if re.match(r'^-?\d+(\.\d*)?$', p)]
                    if len(nums) >= 3:
                        coords.append([float(nums[-3]), float(nums[-2]), float(nums[-1])])
    except FileNotFoundError:
        print(f"ERROR: file not found: {path}")
        return np.zeros((0, 3))

    return np.array(coords)


# ============================================================
# RMSD metrics
# ============================================================
def compute_metrics(A, B):
    """Return RMSD, mean abs error, max deviation."""
    n = min(len(A), len(B))
    if n == 0:
        return None

    A = A[:n]
    B = B[:n]

    diff = A - B
    dists = np.linalg.norm(diff, axis=1)

    return {
        "N compared": n,
        "RMSD": float(np.sqrt((dists ** 2).mean())),
        "Mean abs error": float(dists.mean()),
        "Max deviation": float(dists.max()),
    }


# ============================================================
# MAIN
# ============================================================
input_path = "input.pdb"
full_path  = "output/full_positions.pdb"
cell_path  = "output/cell_positions.pdb"
nbl_path   = "output/nbl_positions.pdb"

files = {
    "input": load_pdb(input_path),
    "full":  load_pdb(full_path),
    "cell":  load_pdb(cell_path),
    "nbl":   load_pdb(nbl_path),
}

# ------------------------------------------------------------
# Atom counts
# ------------------------------------------------------------
print("=== Atom Counts ===")
for name, arr in files.items():
    print(f"{name:6s}: {len(arr)} atoms")
print()

# ------------------------------------------------------------
# Pairwise RMSD comparisons
# ------------------------------------------------------------
pairs = [
    ("input", "full"),
    ("input", "cell"),
    ("input", "nbl"),
    ("full",  "cell"),
    ("full",  "nbl"),
    ("cell",  "nbl"),
]

print("=== Pairwise Differences (RMSD) ===")
for A, B in pairs:
    print(f"\nComparing {A} vs {B}...")
    metrics = compute_metrics(files[A], files[B])
    if metrics is None:
        print("  ERROR: No overlapping atoms!")
        continue

    for k, v in metrics.items():
        print(f"  {k:16s}: {v}")

print("\nDone.\n")
