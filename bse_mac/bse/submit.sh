#!/bin/bash
#SBATCH -A m4868
#SBATCH -C cpu
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time 00:02:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel_chabeda@berkeley.edu

export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

module load PrgEnv-intel/8.5.0
module load cray-fftw/3.3.10.6

#srun -n 16 --cpu-bind=cores -c 16 valgrind --log-file=valgrind-%q{SLURM_PROCID}.dat ./Filter_mpi.x 2> error_val.dat > run_val.dat
srun -n 2 --cpu-bind=cores -c 8 ./BSE.x 2> error_mpi.dat > run_mpi.dat
