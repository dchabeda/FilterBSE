#!/usr/bin/env python3
"""
Passivate perovskite nanocrystals.

Input
-----
conf.par
Format:
N
Sym x y z

Coordinates are assumed to be in BOHR.

Outputs
-------
1. passivated_conf.par   (BOHR)
2. passivated_conf.xyz   (ANGSTROM)

Algorithm
---------
- Detect surface atoms from bounding box proximity
- Determine face / edge / corner atoms
- Determine ligand directions:
      * Halides: opposite Pb-X bond
      * Cs: outward surface normal
      * Pb: opposite average Pb-X directions
- Add:
      face   -> 1 ligand
      edge   -> 2 ligands
      corner -> 3 ligands
- Ligand distance = 0.25 * typical bond length
"""

import numpy as np
from collections import defaultdict


# ============================================================
# CONSTANTS
# ============================================================

BOHR_TO_ANG = 0.52917721092

SURFACE_TOL = 2.0      # bohr
NN_CUTOFF = 8.0        # bohr
LIGAND_SCALE = 0.25

HALIDES = {"I", "Br", "Cl", "F"}
ATOM_TO_LIG = {"I": "PA1", "I0": "PA1", "I1": "PA1", "Cs": "PA2"}
LIG_TO_XYZ = {"PA1": "Be", "PA2": "Li"}

# Typical bond lengths (BOHR)
BOND_LENGTHS = {
    ("Pb", "I"): 6.2,
    ("Pb", "Br"): 5.8,
    ("Pb", "Cl"): 5.4,
    ("Pb", "F"): 4.8,
    ("Cs", "I"): 8.0,
}


# ============================================================
# IO
# ============================================================

def read_confpar(filename):

    atoms = []

    with open(filename, "r") as f:
        lines = [x.strip() for x in f if x.strip()]

    natoms = int(lines[0])

    for line in lines[1:natoms+1]:

        toks = line.split()

        atoms.append({
            "symbol": toks[0],
            "pos": np.array(list(map(float, toks[1:4])))
        })

    return atoms


def write_confpar(filename, atoms):

    with open(filename, "w") as f:

        f.write(f"{len(atoms)}\n")

        for atom in atoms:

            x, y, z = atom["pos"]

            f.write(
                f"{atom['symbol']:4s}"
                f"{x:16.8f}"
                f"{y:16.8f}"
                f"{z:16.8f}\n"
            )


def write_xyz(filename, atoms):

    with open(filename, "w") as f:

        f.write(f"{len(atoms)}\n")
        f.write("passivated nanocrystal\n")

        for atom in atoms:

            x, y, z = atom["pos"] * BOHR_TO_ANG
            if atom["symbol"] in ATOM_TO_LIG.values():
                atom["symbol"] = LIG_TO_XYZ[atom["symbol"]]

            f.write(
                f"{atom['symbol']:4s} "
                f"{x:12.6f} "
                f"{y:12.6f} "
                f"{z:12.6f}\n"
            )


# ============================================================
# UTILITIES
# ============================================================

def unit(v):

    n = np.linalg.norm(v)

    if n < 1e-12:
        return np.zeros(3)

    return v / n


def build_neighbor_list(atoms, cutoff):

    neighs = defaultdict(list)

    positions = np.array([a["pos"] for a in atoms])

    for i in range(len(atoms)):

        for j in range(i+1, len(atoms)):

            rij = positions[j] - positions[i]

            dist = np.linalg.norm(rij)

            if dist < cutoff:

                neighs[i].append(j)
                neighs[j].append(i)

    return neighs


def get_typical_bond_length(sym1, sym2):

    if (sym1, sym2) in BOND_LENGTHS:
        return BOND_LENGTHS[(sym1, sym2)]

    if (sym2, sym1) in BOND_LENGTHS:
        return BOND_LENGTHS[(sym2, sym1)]

    return 6.0


# ============================================================
# SURFACE DETECTION
# ============================================================

def detect_surface_info(atoms):

    positions = np.array([a["pos"] for a in atoms])

    mins = positions.min(axis=0)
    maxs = positions.max(axis=0)

    surface_info = []

    for atom in atoms:

        x, y, z = atom["pos"]

        dirs = []

        if abs(x - mins[0]) < SURFACE_TOL:
            dirs.append(np.array([-1., 0., 0.]))

        if abs(x - maxs[0]) < SURFACE_TOL:
            dirs.append(np.array([1., 0., 0.]))

        if abs(y - mins[1]) < SURFACE_TOL:
            dirs.append(np.array([0., -1., 0.]))

        if abs(y - maxs[1]) < SURFACE_TOL:
            dirs.append(np.array([0., 1., 0.]))

        if abs(z - mins[2]) < SURFACE_TOL:
            dirs.append(np.array([0., 0., -1.]))

        if abs(z - maxs[2]) < SURFACE_TOL:
            dirs.append(np.array([0., 0., 1.]))

        surface_info.append(dirs)

    return surface_info


# ============================================================
# LIGAND DIRECTIONS
# ============================================================

def ligand_directions(i, atoms, neighs, surface_dirs):

    atom = atoms[i]

    sym = atom["symbol"]
    pos = atom["pos"]

    # --------------------------------------------------------
    # HALIDES
    # --------------------------------------------------------

    if sym in HALIDES:

        pb_neighbors = []

        for j in neighs[i]:

            if atoms[j]["symbol"] == "Pb":
                pb_neighbors.append(j)

        if len(pb_neighbors) == 0:
            return surface_dirs

        vec = np.zeros(3)

        for j in pb_neighbors:

            pbpos = atoms[j]["pos"]

            vec += unit(pos - pbpos)

        return [unit(vec)]

    # --------------------------------------------------------
    # Cs
    # --------------------------------------------------------

    elif sym == "Cs":

        return [unit(d) for d in surface_dirs]

    # --------------------------------------------------------
    # Pb
    # --------------------------------------------------------

    elif sym == "Pb":

        vec = np.zeros(3)

        for j in neighs[i]:

            nsym = atoms[j]["symbol"]

            if nsym in HALIDES:

                npos = atoms[j]["pos"]

                vec += unit(pos - npos)

        if np.linalg.norm(vec) < 1e-12:
            return surface_dirs

        return [unit(vec)]

    # --------------------------------------------------------
    # fallback
    # --------------------------------------------------------

    return [unit(d) for d in surface_dirs]


# ============================================================
# PASSIVATION
# ============================================================

def passivate(atoms):

    neighs = build_neighbor_list(atoms, NN_CUTOFF)

    surface_info = detect_surface_info(atoms)

    new_atoms = list(atoms)

    n_added = 0

    for i, dirs in enumerate(surface_info):

        if len(dirs) == 0:
            continue

        atom = atoms[i]

        sym = atom["symbol"]
        pos = atom["pos"]

        # face / edge / corner
        nlig = min(len(dirs), 3)

        lig_dirs = ligand_directions(i, atoms, neighs, dirs)

        # If only one main direction exists,
        # perturb for edges/corners
        if len(lig_dirs) == 1 and nlig > 1:

            base = lig_dirs[0]

            new_dirs = []

            for d in dirs:
                new_dirs.append(unit(base + 0.4*d))

            lig_dirs = new_dirs[:nlig]

        # Typical bond length
        if sym in HALIDES:

            bond_len = get_typical_bond_length("Pb", sym)

        elif sym == "Cs":

            bond_len = get_typical_bond_length("Cs", "I")

        else:

            bond_len = 6.0

        lig_dist = LIGAND_SCALE * bond_len

        for d in lig_dirs[:nlig]:

            lig_pos = pos + lig_dist * unit(d)

            new_atoms.append({
                "symbol": ATOM_TO_LIG[sym],
                "pos": lig_pos
            })

            n_added += 1

    print(f"Added {n_added} ligands")

    return new_atoms


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    atoms = read_confpar("conf.par")

    print(f"Read {len(atoms)} atoms")

    passivated = passivate(atoms)

    write_confpar("passivated_conf.par", passivated)

    write_xyz("passivated_conf.xyz", passivated)

    print("Wrote:")
    print("  passivated_conf.par   (BOHR)")
    print("  passivated_conf.xyz   (ANGSTROM)")
