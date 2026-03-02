import matplotlib.pyplot as plt
import numpy as np
import os

# =============================
# Ensure output folder exists
# =============================
os.makedirs("output", exist_ok=True)

# =============================
# Prepare summary file
# =============================
summary_path = "output/energy_summary.txt"
summary = open(summary_path, "w")
summary.write("=== ENERGY SUMMARY REPORT ===\n\n")

# =============================
# Load energy CSV
# =============================
def load_energy_csv(path):
    if not os.path.exists(path):
        return None

    steps, K, U, E = [], [], [], []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("step"):
                continue

            parts = line.split(",")
            if len(parts) < 4:
                continue

            try:
                steps.append(int(parts[0]))
                K.append(float(parts[1]))
                U.append(float(parts[2]))
                E.append(float(parts[3]))
            except:
                continue

    if len(steps) == 0:
        return None

    return np.array(steps), np.array(K), np.array(U), np.array(E)


# =============================
# File definitions (UPDATED)
# =============================
files = {
    "full":         "../output/full_energies.csv",
    "cell":         "../output/cell_energies.csv",
    "cell_half_mt": "../output/cell_half_mt_energies.csv",
    "nbl":          "../output/nbl_energies.csv",
    "final":        "../output/final_energies.csv"
}

# =============================
# Create main figure (DYNAMIC ROWS)
# =============================
fig, axes = plt.subplots(len(files), 1, figsize=(12, 4 * len(files)), sharex=True)

final_data = None

for ax, (name, path) in zip(axes, files.items()):
    data = load_energy_csv(path)
    if data is None:
        ax.text(0.5, 0.5, f"No data for {name}", ha='center', va='center')
        summary.write(f"{name.upper()}: MISSING FILE ({path})\n\n")
        continue

    steps, K, U, E = data

    if name == "final":
        final_data = data

    # =============================
    # Energy conservation check
    # =============================
    E0 = E[0]
    Ef = E[-1]
    dE = Ef - E0
    pct = (dE / E0) * 100

    # Write to summary
    summary.write(f"--- {name.upper()} ENERGY CHECK ---\n")
    summary.write(f"Initial E0: {E0:.6e}\n")
    summary.write(f"Final Ef : {Ef:.6e}\n")
    summary.write(f"ΔE       : {dE:.6e}\n")
    summary.write(f"ΔE/E0    : {pct:.4f} %\n\n")

    # =============================
    # Plot energies
    # =============================
    ax.plot(steps, K, label="Kinetic", alpha=0.8)
    ax.plot(steps, U, label="Potential", alpha=0.8)
    ax.plot(steps, E, label="Total", linewidth=2.0)

    ax.text(
        0.02, 0.95,
        f"Drift = {pct:.3f}%",
        transform=ax.transAxes,
        fontsize=11,
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )

    ax.set_title(f"{name.upper()} Energy vs Step")
    ax.set_ylabel("Energy")
    ax.grid(True)
    ax.legend()

axes[-1].set_xlabel("Step")

plt.tight_layout()
plt.savefig("output/energy_all_comparison.png", dpi=300)

# =============================
# LOG-PLOT for FINAL
# =============================
if final_data is not None:
    steps, K, U, E = final_data
    eps = 1e-12

    Klog = np.log10(np.abs(K) + eps)
    Ulog = np.log10(np.abs(U) + eps)
    Elog = np.log10(np.abs(E) + eps)

    # summary for log plot
    summary.write("=== FINAL RUN LOG-SCALE ENERGY INFO ===\n")
    summary.write(f"log10|K| range: {np.min(Klog):.4f} to {np.max(Klog):.4f}\n")
    summary.write(f"log10|U| range: {np.min(Ulog):.4f} to {np.max(Ulog):.4f}\n")
    summary.write(f"log10|E| range: {np.min(Elog):.4f} to {np.max(Elog):.4f}\n")
    summary.write(f"K spike magnitude (orders): {np.max(Klog) - np.min(Klog):.4f}\n")
    summary.write(f"U spike magnitude (orders): {np.max(Ulog) - np.min(Ulog):.4f}\n")
    summary.write(f"E flatness (orders): {np.max(Elog) - np.min(Elog):.4f}\n\n")

    plt.figure(figsize=(12, 6))
    plt.plot(steps, Klog, label="log10|K|", alpha=0.8)
    plt.plot(steps, Ulog, label="log10|U|", alpha=0.8)
    plt.plot(steps, Elog, label="log10|E|", linewidth=2.0)

    plt.title("FINAL Simulation: Log-Scale Energies (log10|E|)")
    plt.xlabel("Step")
    plt.ylabel("log10(|Energy|)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("output/energy_final_logplot.png", dpi=300)

summary.close()
print(f"\nEnergy summary written to: {summary_path}")
print("Comparison saved: output/energy_all_comparison.png")
print("Final log-plot saved: output/energy_final_logplot.png\n")