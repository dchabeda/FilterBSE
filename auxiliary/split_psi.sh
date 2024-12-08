# Define N_1 and N_2
N_1=$1 # Number of bytes for psi-holes.dat
N_2=$2 # Number of bytes for psi-elecs.dat

# Split the first N_1 bytes into psi-holes.dat
dd if=psi.dat of=psi-holes.dat bs=1 count=$N_1

# Split the remaining N_2 bytes into psi-elecs.dat
dd if=psi.dat of=psi-elecs.dat bs=1 skip=$N_1 count=$N_2
