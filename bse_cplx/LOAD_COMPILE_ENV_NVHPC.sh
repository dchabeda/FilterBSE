#!/bin/sh
# ============================================================================
#  Environment for the NVHPC + cray-libsci build of bse_cplx WITH GPU offload
#  (produces bse_cplx_gpu.x via the default Makefile). Use this same module set
#  at BUILD time and at RUN time on Perlmutter GPU nodes.
#
#  IMPORTANT: `source` this in your current shell --
#       source LOAD_COMPILE_ENV_NVHPC.sh
#  Do NOT pipe or redirect its output (e.g. `... | tee`); a pipe runs the
#  module commands in a subshell and the PrgEnv swap is silently lost, leaving
#  `cc` as Intel icx (which rejects -mp=gpu). After sourcing, confirm with:
#       echo $PE_ENV        # -> NVIDIA
#       cc --version        # -> nvc / NVIDIA
# ============================================================================

module load PrgEnv-nvidia    # NVHPC compilers behind the Cray `cc` wrapper
module load cray-fftw        # FFTW for the Hartree potential
module load cray-libsci      # BLAS/LAPACK (LAPACKE_zheev), links -lsci_nvidia

# The Cray cc wrapper links cray-libsci and supplies the MPI include paths that
# bare `nvc` lacks -- always build through `cc`, not `nvc` directly.
