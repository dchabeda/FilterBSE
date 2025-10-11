import numpy as np

temp = 300

global AUTOS 
global prefactor
global kB
AUTOS = 2.4188843265857e-17
c = 137.03604
eps0 = (1/(4*np.pi))

prefactor = ( 3*np.pi) * eps0 * c**3 

kB = 3.167e-6
kBT = kB*temp
beta = 1/kBT

def calc_radiative(w, mu2):
  # First we compute all the lifetimes in ns 
  rates = (w**3 * mu2)/prefactor / AUTOS / 1e9
  #print(f"rates = {rates}")
  # Then we compute the radiative rates.
  zero_mask = (rates != 0.0).astype(bool)
  rates = rates[zero_mask]
  # The lifetimes have an inverse dependence on w and |mu|
  # while rates have a polynomial dependence; we average over rates
  
  # Compute the partition function
  boltz_factors = np.exp(- (w) / kBT)
  Q = np.sum(boltz_factors)
  # Compute the thermal average rates and lifetimes
  #print(f"factors times rates: {boltz_factors * rates}")
  avg_rate = (1/Q) * np.sum(boltz_factors[zero_mask] * rates)  
  avg_lifetime = (1/avg_rate) 

  return avg_rate, avg_lifetime


# Get the dipole moment values for each state and their corresponding energy
os_file = np.loadtxt('OS.dat', dtype=np.float64)
state = os_file[:,0]
ene = os_file[:,2]
os2 = np.sum(os_file[:,-6:]**2, axis=1) # Grab sqrt(mu^2) and square it

ene = ene[:100]
os2 = os2[:100]

rate, lifetime = calc_radiative(ene, os2)

print(f"Rate = {rate*1e3:.3f} us^-1\nLifetime = {lifetime:.3f} ns\n")

outfile = open("radiative.dat", "w")
outfile.write(f"Rate = {rate*1e3:.3f} us^-1\nLifetime = {lifetime:.3f} ns\n")
for i in range(os2.shape[0]):
  outfile.write(f"{os2[i]:.4f}\n")
outfile.close()

# Find first large os^2 state
large_os2 = os2[os2 > 40][0]
large_ene = ene[os2 > 40][0]
large_rate = (large_ene**3 * large_os2) / prefactor / AUTOS / 1e9
large_lifetime = 1 / large_rate
print(f"{large_rate * 1e3:.3f}")
print(f"{large_lifetime:.3f}")