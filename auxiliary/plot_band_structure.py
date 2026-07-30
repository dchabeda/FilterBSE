#!/usr/bin/env python3
"""
Plot a band structure from Filter's per-k eval-<k>.dat files.

Only a subset of the states in each eval file are true eigenvalues; the rest are
poorly-converged ghosts. We therefore keep only states whose eigenvalue variance
(sigma_E, column 3) is below a threshold, then build the hole (valence) and
electron (conduction) manifolds outward from the band edges that straddle the
Fermi energy:

    valence band  n  =  the (n+1)-th converged state at or below E_fermi, going DOWN
    conduction band n =  the (n+1)-th converged state above E_fermi,    going UP

Band n is then drawn by connecting that same index-from-the-edge across k-points,
so the band-edge state (n=0) is the VBM / CBM at every k.

The |k| for each k-point is read from a Filter run log ("These are the k_vecs:"
block); if no log is found it is reconstructed from the unit cell + k-grid in
periodic_input.par. A PDF is written to the current directory.

The left y-axis is in eV; a locked twin axis on the right shows the same range in
Hartree. --emin/--emax cap the energy window (in eV, left-axis units).

Usage:
    python plot_band_structure.py [--dir DIR] [--log FILE] [--threshold 1e-4]
                                  [--fermi E] [--nbands N] [--emin EV] [--emax EV]
                                  [--shift] [--out NAME.pdf]
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless / no display
import matplotlib.pyplot as plt

HARTREE_TO_EV = 27.211386245988


# --------------------------------------------------------------------------- #
#  Parsing helpers
# --------------------------------------------------------------------------- #
def parse_eval_file(path):
    """Return list of (energy, sigma) from an eval-<k>.dat file (idx E sigma)."""
    states = []
    with open(path) as f:
        for line in f:
            tok = line.split()
            if len(tok) < 3:
                continue
            try:
                states.append((float(tok[1]), float(tok[2])))
            except ValueError:
                continue
    return states


def find_eval_files(directory):
    """Return {k_index: path} for every eval-<k>.dat in directory."""
    out = {}
    for p in glob.glob(os.path.join(directory, "eval-*.dat")):
        m = re.search(r"eval-(\d+)\.dat$", os.path.basename(p))
        if m:
            out[int(m.group(1))] = p
    return out


def parse_kvecs_from_log(directory, log=None):
    """
    Scan run logs / output for the 'These are the k_vecs:' block and return
    ({k_index: |k|}, source_path). The block looks like:
        These are the k_vecs:
          0.000000 0.000000 0.000000 mag = 0.000000
          0.000000 0.000000 0.019734 mag = 0.019734
          ...
    If `log` is given, only that file is parsed; otherwise the newest candidate
    run log in `directory` that contains the block wins.
    """
    line_re = re.compile(
        r"^\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+mag\s*=\s*([-\d.eE+]+)"
    )
    if log is not None:
        candidates = [log if os.path.isabs(log) else os.path.join(directory, log)]
    else:
        candidates = []
        for pat in ("run*.dat", "*.out", "*.log", "output.dat", "*.stdout"):
            candidates += glob.glob(os.path.join(directory, pat))
        # newest first, so the most recent run wins
        candidates = sorted(set(candidates), key=os.path.getmtime, reverse=True)

    for path in candidates:
        try:
            with open(path, errors="ignore") as f:
                text = f.read()
        except OSError:
            continue
        idx = text.rfind("These are the k_vecs:")
        if idx < 0:
            continue
        kmags = []
        for line in text[idx:].splitlines()[1:]:
            m = line_re.match(line)
            if not m:
                if kmags:        # block ended
                    break
                continue         # tolerate blank lines before the block body
            kmags.append(float(m.group(4)))
        if kmags:
            return {i: km for i, km in enumerate(kmags)}, path
    return None, None


def kvecs_from_cell(directory, nk_expected):
    """
    Fallback: reconstruct |k| for a Monkhorst-style grid along the lattice axes
    from periodic_input.par (a, b, c in Bohr; nk1, nk2, nk3). Mirrors the simple
    1D-style sampling k_i = i * (2 pi / L) used by init_periodic for nk>1 on one
    axis. Only the single varying axis is supported here (the common nanowire /
    1D-path case); otherwise returns None.
    """
    pp = os.path.join(directory, "periodic_input.par")
    if not os.path.exists(pp):
        return None
    vals = {}
    with open(pp) as f:
        for line in f:
            m = re.match(r"\s*(\w+)\s*=\s*([-\d.eE+]+)", line)
            if m:
                vals[m.group(1)] = m.group(2)
    try:
        abc = {"nk1": ("a", int(float(vals["nk1"]))),
               "nk2": ("b", int(float(vals["nk2"]))),
               "nk3": ("c", int(float(vals["nk3"])))}
    except (KeyError, ValueError):
        return None
    # find the single axis being sampled
    varying = [(axis, nk) for (_, (axis, nk)) in abc.items() if nk > 1]
    if len(varying) != 1:
        return None
    axis, nk = varying[0]
    try:
        L = float(vals[axis])
    except (KeyError, ValueError):
        return None
    dk = 2.0 * np.pi / L
    return {i: i * dk for i in range(nk)}


def parse_fermi(directory):
    """Read 'fermiEnergy = X' (Hartree) from input.par; None if absent."""
    ip = os.path.join(directory, "input.par")
    if not os.path.exists(ip):
        return None
    with open(ip) as f:
        for line in f:
            m = re.match(r"\s*fermiEnergy\s*=\s*([-\d.eE+]+)", line)
            if m:
                return float(m.group(1))
    return None


# --------------------------------------------------------------------------- #
#  Band assembly
# --------------------------------------------------------------------------- #
def build_manifolds(eval_files, kmags, e_fermi, threshold):
    """
    Returns (valence, conduction) where each is a dict
        band_index -> list of (kmag, energy)
    built outward from the Fermi level using only states with sigma < threshold.
    """
    valence, conduction = {}, {}
    n_used = 0
    for k_idx, path in sorted(eval_files.items()):
        if k_idx not in kmags:
            continue
        km = kmags[k_idx]
        conv = [E for (E, sig) in parse_eval_file(path) if sig < threshold]
        n_used += len(conv)
        val = sorted([E for E in conv if E <= e_fermi], reverse=True)  # VBM first
        con = sorted([E for E in conv if E > e_fermi])                 # CBM first
        for b, E in enumerate(val):
            valence.setdefault(b, []).append((km, E))
        for b, E in enumerate(con):
            conduction.setdefault(b, []).append((km, E))
    return valence, conduction, n_used


# --------------------------------------------------------------------------- #
#  Plotting
# --------------------------------------------------------------------------- #
def plot(valence, conduction, e_fermi, args, out_path):
    """Left y-axis in eV (the primary axis), right twin y-axis in Hartree.
    Energies are shifted so E_F = 0 when --shift is set. The energy window
    (--emin/--emax, in the left-axis eV units) caps the view."""
    fig, ax = plt.subplots(figsize=(6.0, 5.0))

    # data are in Hartree; the left axis is eV
    e0 = e_fermi if args.shift else 0.0

    def draw(bands, color, label):
        first = True
        for b in sorted(bands):
            if args.nbands is not None and b >= args.nbands:
                continue
            pts = sorted(bands[b])
            x = [p[0] for p in pts]
            y = [(p[1] - e0) * HARTREE_TO_EV for p in pts]   # -> eV
            ax.plot(x, y, "-o", color=color, ms=3, lw=1.0,
                    label=(label if first else None))
            first = False

    draw(valence, "tab:blue", "valence (holes)")
    draw(conduction, "tab:red", "conduction (electrons)")

    # Fermi line (0 if shifted, else E_F in eV)
    ef_y = 0.0 if args.shift else e_fermi * HARTREE_TO_EV
    ax.axhline(ef_y, ls="--", lw=0.8, color="0.4")
    ax.annotate("$E_F$", xy=(ax.get_xlim()[1], ef_y), xytext=(-2, 2),
                textcoords="offset points", ha="right", va="bottom",
                fontsize=9, color="0.4")

    # Energy window (eV, left-axis units). Apply before reading limits for the
    # Hartree twin so both axes share exactly the same physical range.
    if args.emin is not None or args.emax is not None:
        lo, hi = ax.get_ylim()
        ax.set_ylim(args.emin if args.emin is not None else lo,
                    args.emax if args.emax is not None else hi)

    ax.set_xlabel(r"$|\mathbf{k}|$  (Bohr$^{-1}$)")
    ax.set_ylabel(r"$E - E_F$ (eV)" if args.shift else "Energy (eV)")
    ax.set_title(r"Filter band structure ($\sigma_E < %g$)" % args.threshold)
    ax.legend(loc="best", fontsize=9, frameon=False)

    # Right twin axis in Hartree, locked to the eV axis (pure linear rescale).
    ax_ha = ax.twinx()
    lo_ev, hi_ev = ax.get_ylim()
    ax_ha.set_ylim(lo_ev / HARTREE_TO_EV, hi_ev / HARTREE_TO_EV)
    ax_ha.set_ylabel(r"$E - E_F$ (Ha)" if args.shift else "Energy (Ha)")

    fig.tight_layout()
    fig.savefig(out_path)
    print(f"Wrote {out_path}")


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=".", help="directory with eval-*.dat (default: .)")
    ap.add_argument("--log", default=None,
                    help="Filter run log to read k_vecs from (default: newest run*.dat "
                         "in --dir, e.g. run_periodic_nk8_NF_32.dat)")
    ap.add_argument("--threshold", type=float, default=1e-4,
                    help="max sigma_E to accept a state (default: 1e-4)")
    ap.add_argument("--fermi", type=float, default=None,
                    help="Fermi energy in Hartree (default: read input.par)")
    ap.add_argument("--nbands", type=int, default=None,
                    help="max bands per manifold to draw (default: all converged)")
    ap.add_argument("--emin", type=float, default=None,
                    help="lower energy-window limit (eV, left-axis units; "
                         "relative to E_F if --shift)")
    ap.add_argument("--emax", type=float, default=None,
                    help="upper energy-window limit (eV, left-axis units; "
                         "relative to E_F if --shift)")
    ap.add_argument("--shift", action="store_true",
                    help="shift energies so E_F = 0")
    ap.add_argument("--out", default="band_structure.pdf", help="output PDF name")
    args = ap.parse_args()

    directory = args.dir
    eval_files = find_eval_files(directory)
    if not eval_files:
        sys.exit(f"No eval-*.dat files found in {directory!r}")

    kmags, src = parse_kvecs_from_log(directory, args.log)
    if kmags is None:
        kmags = kvecs_from_cell(directory, len(eval_files))
        if kmags is None:
            sys.exit("Could not find k_vecs in any run log, nor reconstruct them "
                     "from periodic_input.par. Provide a Filter log in --dir.")
        print("k_vecs reconstructed from periodic_input.par unit cell.")
    else:
        print(f"k_vecs parsed from {os.path.basename(src)}")

    e_fermi = args.fermi if args.fermi is not None else parse_fermi(directory)
    if e_fermi is None:
        sys.exit("No Fermi energy: pass --fermi or provide input.par with fermiEnergy.")
    print(f"E_Fermi = {e_fermi:.6g} Ha;  sigma threshold = {args.threshold:g}")

    valence, conduction, n_used = build_manifolds(eval_files, kmags, e_fermi,
                                                  args.threshold)
    n_v = sum(len(v) for v in valence.values())
    n_c = sum(len(c) for c in conduction.values())
    print(f"Kept {n_used} converged states across {len(eval_files)} k-points "
          f"({n_v} valence pts, {n_c} conduction pts).")
    if n_used == 0:
        sys.exit("No states passed the variance threshold; raise --threshold.")

    plot(valence, conduction, e_fermi, args, os.path.join(os.getcwd(), args.out))


if __name__ == "__main__":
    main()
