#!/usr/bin/env python3
# plot_dist_potential.py  (no pandas/numpy required)
#
# Examples:
#  1) Pair scatter U(r) at step 1000:
#     python3 plot_dist_potential.py forces.csv --mode pairs --step 1000 --L 8 --rc 2.5
#  2) Radial distribution g(r) at step 1000:
#     python3 plot_dist_potential.py forces.csv --mode gr --step 1000 --L 8 --rc 2.5 --dr 0.02
#  3) Potential energy vs step (from energies.csv):
#     python3 plot_dist_potential.py energies.csv --mode timeU
#  4) Per-particle U at a snapshot:
#     python3 plot_dist_potential.py forces.csv --mode particleU --step 1000 --L 8 --rc 2.5
#
# Notes:
#  - For 'pairs', 'gr', and 'particleU' we need box length L and cutoff rc.
#  - We compute LJ with **force+potential shift**: U_sf = U - Uc - (r-rc)*Fc_over_r.

import argparse, csv, os, math
import matplotlib.pyplot as plt

def read_forces_snapshot(path, step):
    """Return positions (list of (x,y,z)) for a given step from forces.csv"""
    with open(path, newline="") as f:
        r = csv.reader(f)
        header = next(r, None)
        if not header:
            raise SystemExit("Empty forces CSV")
        idx = {h:i for i,h in enumerate(header)}
        need = ["step","x","y","z"]
        for k in need:
            if k not in idx: raise SystemExit(f"forces.csv missing column: {k}")
        pos=[]
        for row in r:
            try:
                s = int(float(row[idx["step"]]))
            except:
                continue
            if s != step:
                continue
            try:
                x = float(row[idx["x"]]); y=float(row[idx["y"]]); z=float(row[idx["z"]])
            except:
                continue
            pos.append((x,y,z))
    if not pos:
        raise SystemExit(f"No positions found for step={step}")
    return pos

def read_energies_csv(path):
    """Read energies.csv with columns: step,K,U,E (no header)."""
    steps=[]; U=[]; E=[]; K=[]
    with open(path, newline="") as f:
        r = csv.reader(f)
        for row in r:
            if len(row) < 4: continue
            try:
                steps.append(int(float(row[0])))
                K.append(float(row[1])); U.append(float(row[2])); E.append(float(row[3]))
            except:
                continue
    if not steps:
        raise SystemExit("No rows parsed from energies.csv")
    return steps, K, U, E

def minimage(d, L):
    """Minimum-image convention for a single component."""
    if d >  0.5*L: d -= L
    elif d < -0.5*L: d += L
    return d

def lj_shifted_U(r, rc):
    """Lennard-Jones potential with force+potential shift at rc, epsilon=sigma=1."""
    if r >= rc: return 0.0
    inv_r2 = 1.0/(r*r)
    inv_r6 = inv_r2*inv_r2*inv_r2
    U = 4.0*(inv_r6*inv_r6 - inv_r6)
    inv_rc2 = 1.0/(rc*rc)
    inv_rc6 = inv_rc2*inv_rc2*inv_rc2
    Uc = 4.0*(inv_rc6*inv_rc6 - inv_rc6)
    Fc_over_r = 24.0*(2.0*inv_rc6*inv_rc6 - inv_rc6)*inv_rc2   # F(rc)/rc
    Usf = U - Uc - (r - rc)*Fc_over_r
    return Usf

def pair_distances_and_potentials(pos, L, rc):
    """Compute all unique pairs i<j: return lists r_list, U_list (shifted)."""
    n = len(pos)
    rs=[]; Us=[]
    for i in range(n):
        xi, yi, zi = pos[i]
        for j in range(i+1, n):
            xj, yj, zj = pos[j]
            dx = minimage(xi-xj, L)
            dy = minimage(yi-yj, L)
            dz = minimage(zi-zj, L)
            r2 = dx*dx + dy*dy + dz*dz
            r = math.sqrt(r2)
            if r < rc:
                rs.append(r)
                Us.append(lj_shifted_U(r, rc))
    return rs, Us

def plot_pairs(rs, Us, step, outp):
    plt.figure(figsize=(8,6))
    plt.scatter(rs, Us, s=8)
    plt.xlabel("r"); plt.ylabel("U_shifted(r)")
    plt.title(f"Pair U(r) at step {step}")
    plt.grid(True, alpha=0.3); plt.tight_layout()
    plt.savefig(outp, dpi=150)
    print(f"✅ Saved pair U(r) scatter: {outp}")

def plot_gr(rs, L, n, dr, outp):
    # Histogram up to rmax = min(rc, L/2)
    if not rs: 
        raise SystemExit("No pair distances within cutoff to build g(r)")
    rmax = min(max(rs), 0.5*L)
    nbins = max(1, int(rmax/dr))
    bins = [i*dr for i in range(nbins+1)]
    counts = [0]*nbins
    for r in rs:
        b = int(r/dr)
        if 0 <= b < nbins:
            counts[b] += 1
    # Normalize: g(r) = (2 * counts) / (N * rho * 4π r^2 dr)
    rho = n / (L*L*L)
    rcenters=[]; g=[]
    for b in range(nbins):
        rlow = b*dr; rhigh = rlow+dr
        rmid = (rlow+rhigh)/2.0
        shell_vol = 4.0*math.pi*(rmid*rmid)*dr
        # counts are unique pairs (i<j). Multiply by 2 to convert to "neighbors per particle"
        denom = n * rho * shell_vol
        val = (2.0*counts[b])/denom if denom>0 else 0.0
        rcenters.append(rmid); g.append(val)
    plt.figure(figsize=(8,5))
    plt.plot(rcenters, g)
    plt.xlabel("r"); plt.ylabel("g(r)")
    plt.title("Radial Distribution Function (single snapshot)")
    plt.grid(True, alpha=0.3); plt.tight_layout()
    plt.savefig(outp, dpi=150)
    print(f"✅ Saved g(r): {outp}")

def plot_timeU(energies_csv, outp):
    steps, K, U, E = read_energies_csv(energies_csv)
    plt.figure(figsize=(12,6))
    plt.plot(steps, U, label='U (potential)')
    plt.plot(steps, E, label='E (total)', linestyle='--')
    plt.xlabel("Step"); plt.ylabel("Energy")
    plt.title("Potential & Total Energy vs Step")
    plt.legend(); plt.grid(True, alpha=0.3); plt.tight_layout()
    plt.savefig(outp, dpi=150)
    print(f"✅ Saved U/E time plot: {outp}")

def plot_particleU(pos, L, rc, outp):
    """Per-particle potential energy at one snapshot: U_i = 0.5 * sum_j U_ij."""
    n = len(pos)
    Ui = [0.0]*n
    for i in range(n):
        xi, yi, zi = pos[i]
        for j in range(i+1, n):
            xj, yj, zj = pos[j]
            dx = minimage(xi-xj, L)
            dy = minimage(yi-yj, L)
            dz = minimage(zi-zj, L)
            r = math.sqrt(dx*dx+dy*dy+dz*dz)
            if r < rc:
                Uij = lj_shifted_U(r, rc)
                Ui[i] += 0.5*Uij
                Ui[j] += 0.5*Uij
    plt.figure(figsize=(10,5))
    plt.plot(range(n), Ui, marker='o', linestyle='None', markersize=4)
    plt.xlabel("particle id"); plt.ylabel("U_i (shifted)")
    plt.title("Per-particle potential energy (snapshot)")
    plt.grid(True, alpha=0.3); plt.tight_layout()
    plt.savefig(outp, dpi=150)
    print(f"✅ Saved per-particle U: {outp}")

def main():
    ap = argparse.ArgumentParser(description="Distance & potential energy analysis")
    ap.add_argument("csv", help="forces.csv (for pairs/gr/particleU) or energies.csv (for timeU)")
    ap.add_argument("--mode", choices=["pairs","gr","timeU","particleU"], default="pairs")
    ap.add_argument("--step", type=int, default=0, help="Snapshot step for pairs/gr/particleU")
    ap.add_argument("--L", type=float, default=8.0, help="Box length L")
    ap.add_argument("--rc", type=float, default=2.5, help="Cutoff radius rc")
    ap.add_argument("--dr", type=float, default=0.02, help="Bin width for g(r)")
    ap.add_argument("--out", default=None, help="Output PNG path (optional)")
    args = ap.parse_args()

    base = os.path.splitext(os.path.basename(args.csv))[0]

    if args.mode == "timeU":
        out = args.out or f"{base}_U_time.png"
        plot_timeU(args.csv, out)
        return

    # modes using forces.csv positions
    pos = read_forces_snapshot(args.csv, args.step)

    if args.mode == "pairs":
        rs, Us = pair_distances_and_potentials(pos, args.L, args.rc)
        if not rs:
            raise SystemExit("No pairs within cutoff; adjust rc.")
        out = args.out or f"{base}_pairs_step{args.step}.png"
        plot_pairs(rs, Us, args.step, out)

    elif args.mode == "gr":
        rs, _ = pair_distances_and_potentials(pos, args.L, args.rc)
        out = args.out or f"{base}_gr_step{args.step}.png"
        plot_gr(rs, args.L, len(pos), args.dr, out)

    else:  # particleU
        out = args.out or f"{base}_particleU_step{args.step}.png"
        plot_particleU(pos, args.L, args.rc, out)

if __name__ == "__main__":
    main()
