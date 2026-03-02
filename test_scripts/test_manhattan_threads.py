#!/usr/bin/env python3
import argparse
import csv
import os
import re
import shutil
import subprocess
from math import sqrt

PDB_X_COL = (30, 38)  # 0-based slice [30:38]
PDB_Y_COL = (38, 46)
PDB_Z_COL = (46, 54)

def read_pdb_xyz(path: str):
    """
    Reads ATOM/HETATM coordinates from a PDB file.
    Returns a list of (x,y,z) in file order.
    """
    coords = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[PDB_X_COL[0]:PDB_X_COL[1]])
                    y = float(line[PDB_Y_COL[0]:PDB_Y_COL[1]])
                    z = float(line[PDB_Z_COL[0]:PDB_Z_COL[1]])
                except ValueError:
                    continue
                coords.append((x, y, z))
    return coords

def compare_coords(a, b):
    """
    Returns: (rms, max_abs, max_r)
    rms: RMS distance per atom
    max_abs: max absolute coordinate difference across all atoms/axes
    max_r: max Euclidean distance for any atom
    """
    if len(a) != len(b):
        raise ValueError(f"Atom count mismatch: {len(a)} vs {len(b)}")

    n = len(a)
    sum_r2 = 0.0
    max_abs = 0.0
    max_r = 0.0

    for (ax, ay, az), (bx, by, bz) in zip(a, b):
        dx = ax - bx
        dy = ay - by
        dz = az - bz

        max_abs = max(max_abs, abs(dx), abs(dy), abs(dz))

        r2 = dx*dx + dy*dy + dz*dz
        sum_r2 += r2
        max_r = max(max_r, sqrt(r2))

    rms = sqrt(sum_r2 / max(1, n))
    return rms, max_abs, max_r

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def run_cmd(cmd, cwd=None):
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return p.returncode, p.stdout

def extract_time(stdout: str, label: str):
    """
    Extracts time printed like:
    'MANHATTAN CELL-LIST (Pthread) done. Time: 0.399399 seconds'
    """
    # label is something like "MANHATTAN CELL-LIST (Pthread) done. Time:"
    pattern = re.escape(label) + r"\s*([0-9]*\.?[0-9]+)\s*seconds"
    m = re.search(pattern, stdout)
    return float(m.group(1)) if m else None

def main():
    ap = argparse.ArgumentParser(description="Test Manhattan accuracy + timing across pthread counts.")
    ap.add_argument("input_pdb", help="Input PDB, e.g. input/random_particles-1024.pdb")
    ap.add_argument("--threads", default="1,2,4,8,16",
                    help="Comma-separated list of thread counts to test (default: 1,2,4,8,16)")
    ap.add_argument("--exe", default="./md_run", help="Path to md_run (default: ./md_run)")
    ap.add_argument("--outdir", default="output/manhattan_tests", help="Directory to store test artifacts")
    ap.add_argument("--ref", choices=["full", "cell", "nbl", "nbl_pthread", "cell_mt"],
                    default="full",
                    help="Which reference output PDB to compare against (default: full)")
    args = ap.parse_args()

    thread_list = [int(x.strip()) for x in args.threads.split(",") if x.strip()]
    if not thread_list:
        raise SystemExit("No threads provided.")

    ensure_dir(args.outdir)

    # Map reference choice -> expected output filename from your C program
    ref_map = {
        "full": "output/full_positions.pdb",
        "cell": "output/cell_positions.pdb",
        "nbl": "output/nbl_positions.pdb",
        "nbl_pthread": "output/nbl_pthread_positions.pdb",
        "cell_mt": "output/cell_half_mt_positions.pdb",
    }
    ref_path = ref_map[args.ref]

    # Expected Manhattan outputs from your main.c
    manh_pdb = "output/manhattan_pthread_positions.pdb"
    manh_csv = "output/manhattan_pthread_energies.csv"

    summary_rows = []
    print(f"\nRunning Manhattan thread sweep: {thread_list}")
    print(f"Reference: {args.ref} ({ref_path})")
    print(f"Artifacts saved to: {args.outdir}\n")

    for t in thread_list:
        print(f"=== RUN: --threads {t} ===")
        cmd = [args.exe, args.input_pdb, "--threads", str(t)]
        code, out = run_cmd(cmd)

        # Save stdout log for this run
        log_path = os.path.join(args.outdir, f"run_t{t}.log")
        with open(log_path, "w") as f:
            f.write(out)

        if code != 0:
            print(f"Run failed for threads={t} (exit {code}). See {log_path}")
            summary_rows.append({
                "threads": t,
                "manhattan_time_s": "",
                "rms_dist": "",
                "max_abs_coord_diff": "",
                "max_atom_dist": "",
                "status": f"FAIL(exit={code})",
                "log": log_path
            })
            continue

        # Check files exist
        if not os.path.exists(manh_pdb):
            print(f"Missing {manh_pdb}. Manhattan may not be enabled. See {log_path}")
            summary_rows.append({
                "threads": t,
                "manhattan_time_s": "",
                "rms_dist": "",
                "max_abs_coord_diff": "",
                "max_atom_dist": "",
                "status": "FAIL(missing_manhattan_pdb)",
                "log": log_path
            })
            continue

        if not os.path.exists(ref_path):
            print(f"Missing reference PDB {ref_path}. See {log_path}")
            summary_rows.append({
                "threads": t,
                "manhattan_time_s": "",
                "rms_dist": "",
                "max_abs_coord_diff": "",
                "max_atom_dist": "",
                "status": "FAIL(missing_reference_pdb)",
                "log": log_path
            })
            continue

        # Copy outputs so they don't get overwritten by the next run
        manh_pdb_copy = os.path.join(args.outdir, f"manhattan_t{t}.pdb")
        ref_pdb_copy  = os.path.join(args.outdir, f"{args.ref}_t{t}.pdb")
        shutil.copyfile(manh_pdb, manh_pdb_copy)
        shutil.copyfile(ref_path, ref_pdb_copy)

        if os.path.exists(manh_csv):
            shutil.copyfile(manh_csv, os.path.join(args.outdir, f"manhattan_t{t}.csv"))

        # Compute accuracy metrics
        manh_xyz = read_pdb_xyz(manh_pdb_copy)
        ref_xyz  = read_pdb_xyz(ref_pdb_copy)

        try:
            rms, max_abs, max_r = compare_coords(manh_xyz, ref_xyz)
            status_str = "OK"
        except Exception as e:
            rms, max_abs, max_r = "", "", ""
            status_str = f"FAIL(compare_error={e})"

        # Extract Manhattan time from stdout
        manh_time = extract_time(out, "MANHATTAN CELL-LIST (Pthread) done. Time:")

        print(f"  manhattan_time_s: {manh_time}")
        print(f"  rms_dist:         {rms}")
        print(f"  max_abs_diff:     {max_abs}")
        print(f"  max_atom_dist:    {max_r}")
        print(f"  log: {log_path}\n")

        summary_rows.append({
            "threads": t,
            "manhattan_time_s": manh_time if manh_time is not None else "",
            "rms_dist": rms,
            "max_abs_coord_diff": max_abs,
            "max_atom_dist": max_r,
            "status": status_str,
            "log": log_path
        })

    # Write summary CSV
    summary_path = os.path.join(args.outdir, "summary.csv")
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["threads", "manhattan_time_s", "rms_dist", "max_abs_coord_diff", "max_atom_dist", "status", "log"]
        )
        w.writeheader()
        for row in summary_rows:
            w.writerow(row)

    print(f"Done. Summary written to: {summary_path}")

if __name__ == "__main__":
    main()