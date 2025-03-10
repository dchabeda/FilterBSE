import numpy as np

FERMI = -0.19
TOL = 1e-7  # Tolerance to identify degeneracies

# Load the data
data = np.loadtxt("eval.par")
indices, energies, variances = data[:, 0], data[:, 1], data[:, 2]

# Split into hole and electron energies
hole_mask = energies < FERMI
elec_mask = energies > FERMI

hole_energies = data[hole_mask][::-1]  # Sort hole energies in decreasing order
elec_energies = data[elec_mask]

# Function to group degenerate states
def find_degeneracies(states):
    groups = []
    current_group = [states[0]]

    for i in range(1, len(states)):
        if abs(states[i, 1] - current_group[-1][1]) < TOL:
            current_group.append(states[i])
        else:
            groups.append(current_group)
            current_group = [states[i]]
    groups.append(current_group)

    return groups

# Group degenerate states
hole_groups = find_degeneracies(hole_energies)
elec_groups = find_degeneracies(elec_energies)

# Write the degeneracy levels to file
def write_degeneracies(groups, filename):
    with open(filename, "w") as f:
        num_states = 0
        for group in groups:
            num_states += len(group)
            f.write(f"{num_states}\n")

write_degeneracies(hole_groups, "hole_degen.dat")
write_degeneracies(elec_groups, "elec_degen.dat")

print("Degeneracy files written.")
