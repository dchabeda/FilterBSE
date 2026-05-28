import numpy as np
import sys

if len(sys.argv) != 5:
  print("\nusage: python make_gaussian_lig_pot.py [a] [b] [dr] [filename]\n")
  exit()

a = float(sys.argv[1])
b = float(sys.argv[2])
r_max = float(sys.argv[3])
filename = sys.argv[4]

# First and last distance in Bohr
r_initial = 0.0
n_points = 4096

r = np.linspace(r_initial, r_max, n_points)
pot = a * np.exp(-b * r**2)

output = np.concatenate((r.reshape(-1, 1), pot.reshape(-1, 1)), axis=1)

np.savetxt(f"{filename}", output, fmt = '%.8f %.8f')

print(f"Done making {filename} with a = {a}, b = {b}, r_max = {r_max}")
