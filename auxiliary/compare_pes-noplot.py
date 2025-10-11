import numpy as np
# import matplotlib.pyplot as plt
import os
import argparse
#from lammps import lammps

parser = argparse.ArgumentParser()
parser.add_argument('Nmode', type=int)
# parser.add_argument('max_step', type=float)
parser.add_argument('nstep', type=int)
# =============================================================================
# Define related functions...
# =============================================================================

TO_A = np.sqrt(1.60218e-19)*1e-15/np.sqrt(1.66054e-27) * 1e10 # convert to Angstrom
TO_B = TO_A * 1.88973 # convert to Bohr

def read_min_geom_data(filename):
    with open(filename, 'r') as f:
        data = f.read()
    return data


def parse_min_geom_data(data):
  # Get the min geom coords from init-min.xyz
  init_min = np.loadtxt('init-min.xyz', skiprows=2, dtype='object')
  init_min[:, 0] = init_min[:, 0].astype(np.int32)
  init_min[:, 1:] = init_min[:, 1:].astype(np.float64)
  
  # Get the lammps input data from argument [data]
  lines = data.strip().split('\n')
  
  # Get the number of atoms from the correct line
  N = int(lines[2].split()[0])
  
  # Locate the Atoms section
  atom_start_idx = None
  for idx, line in enumerate(lines):
      if "Atoms" in line:
          atom_start_idx = idx + 2  # The atom data starts two lines after "Atoms"
          break
  
  # If we can't find the Atoms section, raise an error
  if atom_start_idx is None:
      raise ValueError("Cannot locate Atoms section in the LAMMPS data.")
  
  # Extract atom data
  atom_data = []
  for i in range(atom_start_idx, atom_start_idx + N):
      if i >= len(lines):  # Check to ensure we're not reading past the end of the file
          raise ValueError("Unexpected end of LAMMPS data while reading atoms.")
      
      tokens = lines[i].split()
      if len(tokens) < 7:  # Ensure each line has the required number of columns
          raise ValueError(f"Unexpected line format in Atoms section: {lines[i]}")
      
      atom_id, mol_id, atom_type, q, x, y, z = int(tokens[0]), int(tokens[1]), int(tokens[2]), float(tokens[3]), float(tokens[4]), float(tokens[5]), float(tokens[6])
      idx = i - atom_start_idx
      if init_min[idx, 0] == atom_type:
        x_min = init_min[idx, 1]
        y_min = init_min[idx, 2]
        z_min = init_min[idx, 3]
      else:
        print(f"Error: mismatch between plmp atom type {atom_type} and init-min {init_min[idx, 0]}")
        exit()

      atom_data.append([atom_id, mol_id, atom_type, q, x_min, y_min, z_min])
  
  return N, np.array(atom_data)

def save_xyz(filename, atom_data):
    with open(filename, 'w') as f:
        f.write(str(len(atom_data)) + '\n')
        f.write('\n')
        for n in range(len(atom_data)):
            if atom_data[n,2]==1:
                atm_type = "Cs"
            elif atom_data[n,2]==2:
                atm_type = "Pb"
            elif atom_data[n,2]==3:
                atm_type = "I"
            elif atom_data[n,2]==4:
                atm_type = "Cs"
            else:
                print("Error in reading atom type...")
                exit()

            f.write(atm_type + ' ' + str(atom_data[n,-3])+ ' ' + str(atom_data[n,-2])+ ' ' + str(atom_data[n,-1]) + '\n')
        f.close()


def displace_atoms(atom_data, normal_mode, disp):
    atom_data_cp = atom_data.copy()
    atom_data_cp[:, 4:] += normal_mode * disp / np.sqrt(mass.reshape(-1,1)) * TO_A
    # print("**********")
    # print(atom_data_cp[:10])
    # print("**********")
    return atom_data_cp

def write_displaced_files(filename, displaced_data):
    with open("plmp.dat", 'r') as f:
        lines = f.readlines()

    # Write to a new file with the same headers and replaced atomic data
    with open("displaced_" + filename, 'w') as f:
        # Write header
        for line in lines[:19]:
            f.write(line)
        # f.write("\n") # delete this if you only have 2 atom type
        # Write displaced atomic data
        for row in displaced_data:
            f.write(f"{int(row[0])}  {int(row[1])}  {int(row[2])}  {row[3]:.8g} {row[4]:.8g} {row[5]:.8g} {row[6]:.8g}\n")


def read_mass(natoms):
    filename = 'conf.par'
    if not os.path.exists(filename):
        print('Cannot find ' + filename + '! Exiting...\n')
        exit()

    lines = np.array([line.strip().split() for line in open(filename, 'r')][1:])
    at_type = np.array([line[0] for line in lines])[:natoms]
    
    # make sure no passivation atoms
    if ('P1' == at_type).any() or ('P2' == at_type).any():
        print('Error reading ' + filename + '! Exiting...\n')
        exit()

    mass = np.empty(natoms, dtype=float)
    cd_idx = np.where(at_type == 'Cd')[0]
    se_idx = np.where(at_type == 'Se')[0]
    s_idx = np.where(at_type == 'S')[0]
    in_idx = np.where(at_type == 'In')[0]
    as_idx = np.where(at_type == 'As')[0]
    p_idx = np.where(at_type == 'P')[0]
    cs_idx = np.where(at_type == 'Cs')[0]
    pb_idx = np.where(at_type == 'Pb')[0]
    i_idx = np.where(at_type == 'I')[0]
    if (cd_idx.shape[0] + se_idx.shape[0] + s_idx.shape[0] + in_idx.shape[0] + as_idx.shape[0] + p_idx.shape[0] + cs_idx.shape[0] + pb_idx.shape[0] + i_idx.shape[0] != natoms):
        print('Error reading ' + filename + '! Exiting...\n')
        exit()
    # masses in gram/mole
    mass[cd_idx] = 112.411
    mass[se_idx] = 78.960
    mass[s_idx] = 32.065
    mass[in_idx] = 114.818
    mass[as_idx] = 74.922
    mass[p_idx] = 30.974
    mass[cs_idx] = 132.9054
    mass[pb_idx] = 207.2
    mass[i_idx] = 126.90447
    # mass *= (1e-3/6.022140857e23) # kg

    return mass



# Constants
to_hz = 98.22704446943315 / (2.*np.pi)
hbar = 0.6582119569  # eV⋅fs
kB = 8.617333262145e-5 # Boltzmann constant (eV/K)


temp = 300.0 # K
kT = kB*temp # eV
eig = np.loadtxt('eig.dat')
omega = np.loadtxt("w.dat")[:,1][6:]*2*np.pi*10**(12)*10**(-15) # unit fs^-1

args = parser.parse_args()
Nmode = args.Nmode
Q_max = np.sqrt(kT)/omega[Nmode] # sqrt(eV)/(fs^-1) = sqrt(eV)*fs

max_step = Q_max # positive value
nstep = int(args.nstep / 2)
dQ = round(max_step / nstep, 2)
#nstep = int(max_step/dQ) 
disp_list = [n * dQ for n in range(nstep+1)]
disp_list = np.array(disp_list[::-1] + disp_list[1:])

# =============================================================================
# Generate displacement configurations for lammps...
# =============================================================================

# Read the configuration data from external file
lmp_conf = read_min_geom_data('plmp.dat')
natoms, atom_data = parse_min_geom_data(lmp_conf)

normal_mode = eig[:,Nmode].reshape(-1,3)
num_modes = eig.shape[1]
mass = read_mass(natoms)
sr_mass3 = np.sqrt(np.repeat(mass, 3)) # sqrt(g/mol)

for n, disp in enumerate(disp_list):
    displaced_data = displace_atoms(atom_data, normal_mode, disp)
    write_displaced_files("lammpsconf_nx_{0}.par".format(n), displaced_data)

    if n==0:
        save_xyz("max_disp.xyz", displaced_data)

# =============================================================================
# Use lammps SW to evaluate the potential energies for each config...
# =============================================================================
#os.system(f"rm pes_output_mode_{Nmode}_dx_{dQ}.dat") # avoid re-appending to the previous file
#os.system("rm displaced_lammpsconf_nx_*.par")

command = ''

for n, disp in enumerate(disp_list):
    
  command += \
  """
  # CsPbI3 Perovskite nanocubes with LJ+coulomb pair potential

  # variables and units
  units       metal
  variable    time_step index 0.001           # 0.001 ps = 1.0 fs
  variable    start_temp index 0.001
  variable    seed index 359592               # random seed
  dimension   3
  boundary    f f f
 	
  # Displacement and phonon mode variables
  variable disp equal "%.4f"
  variable    Nmode index %d

  atom_style  full

  # file with the atom charges, positions, etc. 
  read_data displaced_lammpsconf_nx_%d.par

  # Force field
  pair_style     lj/cut/coul/cut 500.0 500.0

  # 1 - Cs / 2 - Pb / 3 - I / 4 - Cs
  # set epsilon, sigma
  pair_coeff              1 1 0.07728  3.584 # Cs
  pair_coeff              1 2 0.009161 3.397 # Cs-Pb
  pair_coeff              1 3 0.070267 3.799 # Cs-I

  pair_coeff              2 2 0.001086   3.210 # Pb
  pair_coeff              2 3 0.008330   3.612 # Pb-I

  pair_coeff              3 3 0.06389 4.014 # I

  pair_coeff              1 4 0.07728  3.584 # Cs
  pair_coeff              2 4 0.009161 3.397 # Pb-Cs
  pair_coeff              3 4 0.070267 3.799 # I-Cs
  pair_coeff              4 4 0.07728  3.584 # Cs-Cs

  dielectric 1.0

  # set up neighbor list information
  neighbor    2.0  bin
  neigh_modify  one 1000 check yes delay 0

  # set up the thermodynamic printing in the minimization
  thermo 			  5 	# outputs the thermo info every 5 steps
  thermo_style	custom step temp press pe ke etotal

  # write_dump all xyz disp_nmode_${Nmode}_${disp}.xyz
  
  # Compute energy
  run 0 post no
  
  # Output energy for this displacement
  variable my_pe equal pe
  shell echo ${disp} ${my_pe} >> pes_output_mode_${Nmode}_dx_%.2f.dat

  clear
  """%(disp,Nmode,n,dQ)
  
  #print(command)
lmp_input = open(f"disp_nmode_{Nmode}.in", 'w')
lmp_input.write(command)
lmp_input.close()
  
  # lmp_init = lammps()
  # lmp_init.commands_string(command)
  # lmp_init.close()

# data = np.loadtxt("pes_output_mode_%d_dx_%.4f.dat"%(Nmode,dQ))
# Q = data[:, 0] # Ang
# energies = data[:, 1]
# energies = energies - np.min(energies)

# # phonon modes from harmonic app.
# energies_harm = 0.5 * omega[Nmode]**2 * Q**2
# np.savetxt("harm_pes_output_mode_%d_dx_%.4f.dat"%(Nmode,dQ), np.column_stack((Q, energies_harm)))



