#!/usr/bin/env python3
# plot_forces.py  (pandas-free)
#
# Usage:
#   python3 plot_forces.py forces.csv --mode snapshot --step 500 --scale 0.2
#   python3 plot_forces.py forces.csv --mode time
#   python3 plot_forces.py forces.csv --mode particle --pid 42
#
# Expects header: step,id,fx,fy,fz,x,y,z,vx,vy,vz

import argparse, csv, os, math
import matplotlib.pyplot as plt

def load_forces(path):
    with open(path, newline="") as f:
        r = csv.reader(f)
        try:
            header = next(r)
        except StopIteration:
            raise SystemExit("Empty CSV.")
        idx = {h:i for i,h in enumerate(header)}
        required = ["step","id","fx","fy","fz","x","y","z"]
        missing = [k for k in required if k not in idx]
        if missing:
            raise SystemExit(f"CSV missing columns: {missing}")

        rows = []
        for row in r:
            if not row or len(row) < len(idx):
                continue
            try:
                rows.append({
                    "step": int(float(row[idx["step"]])),
                    "id":   int(float(row[idx["id"]])),
                    "fx":   float(row[idx["fx"]]),
                    "fy":   float(row[idx["fy"]]),
                    "fz":   float(row[idx["fz"]]),
                    "x":    float(row[idx["x"]]),
                    "y":    float(row[idx["y"]]),
                    "z":    float(row[idx["z"]]),
                })
            except Exception:
                # skip malformed line
                continue
        if not rows:
            raise SystemExit("No valid data rows parsed from CSV.")
        return rows

def plot_snapshot(rows, step, scale, outpath):
    snap = [r for r in rows if r["step"] == step]
    if not snap:
        raise SystemExit(f"No rows found for step={step}")
    xs  = [r["x"]  for r in snap]
    ys  = [r["y"]  for r in snap]
    fxs = [r["fx"] for r in snap]
    fys = [r["fy"] for r in snap]

    plt.figure(figsize=(10,8))
    plt.quiver(xs, ys, fxs, fys, angles='xy', scale_units='xy',
               scale=(1.0/scale if scale != 0 else 1.0), width=0.002)
    plt.title(f'Force Vectors (XY) at step {step}')
    plt.xlabel('x'); plt.ylabel('y')
    plt.axis('equal'); plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    print(f"✅ Saved snapshot quiver: {outpath}")

def plot_time(rows, outpath):
    by_step = {}
    for r in rows:
        F = math.sqrt(r["fx"]**2 + r["fy"]**2 + r["fz"]**2)
        by_step.setdefault(r["step"], []).append(F)

    steps = sorted(by_step.keys())
    fmean = [sum(by_step[s])/len(by_step[s]) for s in steps]
    fmax  = [max(by_step[s]) for s in steps]

    plt.figure(figsize=(12,6))
    plt.plot(steps, fmean, label='mean |F|')
    plt.plot(steps, fmax,  label='max |F|', linestyle='--')
    plt.title('Force Magnitudes vs Step')
    plt.xlabel('Step'); plt.ylabel('|F|')
    plt.legend(); plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    print(f"✅ Saved time series plot: {outpath}")

def plot_particle(rows, pid, outpath):
    track = [(r["step"], math.sqrt(r["fx"]**2 + r["fy"]**2 + r["fz"]**2))
             for r in rows if r["id"] == pid]
    if not track:
        raise SystemExit(f"No rows for particle id={pid}")
    track.sort()
    steps = [t[0] for t in track]
    mags  = [t[1] for t in track]

    plt.figure(figsize=(12,5))
    plt.plot(steps, mags)
    plt.title(f'|F| over Time for Particle {pid}')
    plt.xlabel('Step'); plt.ylabel('|F|')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    print(f"✅ Saved particle force plot: {outpath}")

def main():
    ap = argparse.ArgumentParser(description="Plot forces from forces.csv")
    ap.add_argument("csv", help="Path to forces.csv")
    ap.add_argument("--mode", choices=["snapshot","time","particle"], default="snapshot")
    ap.add_argument("--step", type=int, default=0, help="Step to plot for snapshot mode")
    ap.add_argument("--scale", type=float, default=0.2, help="Arrow scale for snapshot (bigger=longer)")
    ap.add_argument("--pid", type=int, default=0, help="Particle ID for particle mode")
    ap.add_argument("--out", default=None, help="Output PNG path (optional)")
    args = ap.parse_args()

    rows = load_forces(args.csv)
    base = os.path.splitext(os.path.basename(args.csv))[0]

    if args.mode == "snapshot":
        out = args.out or f"{base}_forces_step{args.step}.png"
        plot_snapshot(rows, args.step, args.scale, out)
    elif args.mode == "time":
        out = args.out or f"{base}_forces_time.png"
        plot_time(rows, out)
    else:
        out = args.out or f"{base}_forces_particle{args.pid}.png"
        plot_particle(rows, args.pid, out)

if __name__ == "__main__":
    main()
