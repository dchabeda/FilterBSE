import os
import glob
import matplotlib.pyplot as plt

# Conversion constants
HARTREE_TO_EV = 27.211386
EV_TO_MEV = 1000.0

output_file = "exciton_convergence.dat"
results = []

for path in glob.glob("bse_*_*"):
    exciton_file = os.path.join(path, "exciton.dat")
    if not os.path.isfile(exciton_file):
        continue

    # Extract n from directory name (bse_n_n)
    try:
        n_str = path.split("_")[1]
        n = int(n_str)
    except Exception:
        continue

    with open(exciton_file, "r") as f:
        lines = f.readlines()

    try:
        e0_hartree = float(lines[1].split()[1])  # row 1 col 1
        e1_hartree = float(lines[2].split()[1])  # row 2 col 1
        ex0_hartree = float(lines[1].split()[4])  # row 1 col 4
        ex1_hartree = float(lines[2].split()[4])  # row 1 col 4
    except (IndexError, ValueError) as e:
        print(f"Error reading {exciton_file}: {e}")
        continue

    # Convert to eV
    e0_ev = e0_hartree * HARTREE_TO_EV
    e1_ev = e1_hartree * HARTREE_TO_EV
    gap_mev = (e1_ev - e0_ev) * EV_TO_MEV

    ex0_ev = ex0_hartree * HARTREE_TO_EV * EV_TO_MEV
    ex1_ev = ex1_hartree * HARTREE_TO_EV * EV_TO_MEV

    results.append((n, e0_ev, e1_ev, gap_mev, ex0_ev, ex1_ev))

# Sort results numerically
results.sort(key=lambda x: x[0])

# Write summary file
with open(output_file, "w") as out:
    out.write("n_elec_hole\tE_0(eV)\tE_1(eV)\tGap(meV)\n")
    for n, e0, e1, gap, ex0, ex1 in results:
        out.write(f"{n}\t{e0:.4f}\t{e1:.4f}\t{gap:.1f}\t{ex0:.2g}\t{ex1:.2g}\n")

print(f"Summary written to {output_file}")

# --- Plot ---
n_vals = [r[0] for r in results]
E0_vals = [r[1] for r in results]
E1_vals = [r[2] for r in results]
Gap_vals = [r[3] for r in results]
Ex0_vals = [r[4] for r in results]
Ex1_vals = [r[5] for r in results]

fig, ax1 = plt.subplots(figsize=(6,4))

# Plot energies
ax1.plot(n_vals, E0_vals, "o-", label="E0 (eV)")
ax1.plot(n_vals, E1_vals, "s-", label="E1 (eV)")
ax1.set_xlabel("n_elec_hole")
ax1.set_ylabel("Energy (eV)")
ax1.legend(loc="upper left")

# Second y-axis for gap
ax2 = ax1.twinx()
ax2.plot(n_vals, Gap_vals, "d--", color="red", label="Gap (meV)")
ax2.set_ylabel("Gap (meV)", color="red")
ax2.tick_params(axis="y", labelcolor="red")

fig.tight_layout()
plt.savefig("exciton_plot.png", dpi=200)
plt.show()

# Plot the exchange energy
fig, ax1 = plt.subplots(figsize=(6,4))

# Plot energies
ax1.plot(n_vals, Ex0_vals, "o-", color="k", label="Ex0")
ax1.set_xlabel("n_elec_hole")
ax1.set_ylabel("Exchange energy (meV)")
ax1.legend(loc="upper left")

# Second y-axis for gap
ax2 = ax1.twinx()
ax2.plot(n_vals, Ex1_vals, "s-", color='r', label="Ex1")
ax2.set_ylabel("Exchange energy (meV)", color="red")
ax2.tick_params(axis="y", labelcolor="red")

fig.tight_layout()
plt.savefig("exchange_plot.png", dpi=200)
plt.show()