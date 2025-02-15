To run the code, first compute filter and obtain output.dat
Create an equilibrium geometry file (same from filter) and label conf_equil.par
Then, create a distorted geometry and label this file conf.par
You will need potCs.par, potPb.par, potI.par, SO_Cs.par, SO_Pb.par, SO_I.par
Modify input.par to meet quasiparticle matrix element needs
Output file pot_coupling_equil.dat has potential matrix elements of the equil potential
File pot_coupling_dist.dat has mat elems of the potential of the distorted geometry
  (The first 4 columns are i energy_i j energy_j, 
  and the last two columns are the re and im part of the integral)
