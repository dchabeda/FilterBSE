# Daniel Chabeda 12.14.25
#

# This directory contains the necessary files to run a
# filter diagonalization job. The configuration is a
# 1 unit cell cube of CsPbI3. Pre-fitted pseudopotentials
# are provided to compute the electronic structure.
# The output of the calculation are quasiparticle 
# eigenstates of the pseudopotential Hamiltonian. These
# states are often utilized as the basis for solving a
# Bethe-Salpether equation for correlated excitons.
# To run BSE -> go to FilterBSE/bse_cplx/EXAMPLES directory

# To run the job, first ensure that you have compiled the
# Filter_mpi.x executable in the parent directory. The
# executable requires Intel MKL, FFTW3, and OpenMP libraries.
# The code is MPI aware and is optimal on MPI-capable systems.

# The ./filter directory contains all of the input files 
# required to run the job EXCEPT for the executable itself,
# Filter_mpi.x. Copy the contents of the ./filter directory
# to a local scratch directory, navigate to the directory on
# a compute-node, and run the executable with:
# srun -n N -c C ./Filter_mpi.x | tee -a run.dat

# If you are not on a SLURM system, you can use
# mpirun -n N ./Filter_mpi.x | tee -a run.dat

# This will run the executable on N MPI ranks and
# stream the output to your screen and a file run.dat
# The number of MPI ranks should be set as a factor of 
# the variable nFilterCycles in input.par. For example,
# with nFilterCycles = 16, fitting values are 16 (optimal), 
# 8, 4, 2, 1. The remaining CPUs available on your machine 
# after designating MPI ranks should be used for multi-threading
# with OpenMP. The Hamiltonian evaluation is parallelized with
# OMP threads, so set the hamThreads parameter to be equal to the
# remaining CPUs-per-rank. If your machine is dual-core, you might
# get better performace by placing only 1 thread on 2 physical CPUs.
# Ex. On a 256 dual-core CPU machine like Perlmutter, assigning 16 MPI
# ranks leaves each rank with 16 CPUs. I can use 16 CPUS / 2 CPUs-per-thread
# = 8 hamThreads: srun -n 16 -c 16 ./Filter_mpi.x

# * Extension: this directory also contains BASH scripts to automate
# the submission of filter jobs. Configure the paths and settings
# to your HPC system and try submitting the job with submit.sh