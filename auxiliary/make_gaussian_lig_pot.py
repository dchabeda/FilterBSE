import numpy as np
import sys

if len(sys.argv) != 5:
  print("\nusage: python make_gaussian_lig_pot.py [a] [b] [dr] [filename]\n")
  exit()

a = float(sys.argv[1])
b = float(sys.argv[2])
dr = float(sys.argv[3])
filename = sys.argv[4]



# First and last distance in Bohr
r_initial = 0.0
r_final = 25.0

r = np.arange(r_initial, r_final, dr)
pot = a * np.exp(-r**2/b)

output = np.concatenate((r.reshape(-1, 1), pot.reshape(-1, 1)), axis=1)

np.savetxt(f"{filename}", output, fmt = '%.8f %.8f')

print(f"Done making {filename} with a = {a}, b = {b}, dr = {dr}")
