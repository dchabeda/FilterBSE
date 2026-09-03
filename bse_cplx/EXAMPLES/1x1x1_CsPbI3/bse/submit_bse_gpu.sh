#!/bin/bash
#SBATCH -J bse_gpu_1x1x1
#SBATCH -o %j-%N.out
#SBATCH -C gpu                    # Perlmutter GPU nodes (4x A100 40GB each)
#SBATCH -q debug                  # debug: <= 30 min; use -q regular for production
#SBATCH -A m4868                  # NERSC allocation (change to your GPU account)
#SBATCH -N 1                      # nodes
#SBATCH --ntasks-per-node=4       # 1 MPI rank per GPU
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=map_gpu:0,1,2,3
#SBATCH -c 32                     # CPU cores per rank
#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel@qernelzoo.com

# --- Runtime environment: the SAME module set used to build bse_cplx_gpu.x ---
source /global/common/software/m4868/FilterBSE/bse_cplx/LOAD_COMPILE_ENV_NVHPC.sh

# One GPU per rank; per-rank device-memory budget for the resident state set.
# BSE_GPU_MEM_GB = A100_GB * fraction / ranks_per_gpu  (40GB A100, 1 rank/GPU -> ~36).
# If the state set exceeds this, that kernel falls back to the CPU automatically.
export BSE_GPU_MEM_GB=36

# Host-side OpenMP (Hartree FFT threads, optical matrix elements, etc.)
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# GPU offload is requested via `gpuAccel = 1` in input.par (the default) and is
# auto-gated by device availability. To force the CPU path for an A/B check,
# either set `gpuAccel = 0` in input.par or uncomment the next line:
# export OMP_TARGET_OFFLOAD=DISABLED

# Optional: verify offload actually lands on the device (verbose).
# export NVCOMPILER_ACC_DEBUG=1

BSE_EXE=/global/common/software/m4868/FilterBSE/bse_cplx/bse_cplx_gpu.x

# Run from this directory (must hold input.par, conf.dat, output.dat).
# NOTE: 1x1x1 is a tiny validation case -- the GPU offload only *wins* at
# production grid/basis sizes; here it exists to check correctness, not speed.
srun ${BSE_EXE} 2> error.dat > run.dat
wait
