# FilterBSE

A suite of electronic-structure codes for semiconductor nanostructures. FilterBSE
computes single-particle (quasiparticle) states by **filter diagonalization** and
optical / excitonic properties by solving the **Bethe–Salpeter equation (BSE)**,
using a semiempirical pseudopotential Hamiltonian on a real-space grid.

📖 **Documentation: https://dchabeda.github.io/FilterBSE/**

## What it does

FilterBSE runs a two-step pipeline:

1. **Filter** — compute quasiparticle states and energies on a real-space grid
   (Chebyshev/Newton filter diagonalization).
2. **BSE** — build the electron–hole kernel (direct + exchange Coulomb) from those
   states and diagonalize it for excitons, oscillator strengths, and optical
   spectra.

The codes are written in C with MPI + OpenMP, and use MKL or Cray libsci
(LAPACK/ScaLAPACK) and FFTW. The BSE solver additionally supports distributed
(ScaLAPACK) eigensolves and optional GPU offload of the Coulomb kernel.

## Repository layout

| Directory | Description |
|---|---|
| `filter_mpi/` | MPI filter-diagonalization (distributed filter → orthogonalize → diagonalize) |
| `filter/` | Serial / reference filter code |
| `periodic_filter/` | Periodic-system filter |
| `bse_cplx/` | Complex Bethe–Salpeter solver (MPI + ScaLAPACK + optional GPU) |
| `bse_real/`, `bse_gpu/` | Real-valued and GPU BSE variants |
| `pots/` | Pseudopotential data |
| `auxiliary/` | Helper tools |
| `docs/` | Documentation source (MkDocs Material) |

The two primary, documented codes are **`filter_mpi`** and **`bse_cplx`**.

## Quick start

Both codes are built with per-directory Makefiles. On Perlmutter, load the build
environment first.

### 1. Filter (`filter_mpi`)

```sh
source LOAD_COMPILE_ENV.sh          # PrgEnv-intel + cray-fftw
cd filter_mpi && make               # produces Filter_mpi.x

# Run in a directory with input.par, conf.par, pot<Sym>.par, SO_<Sym>.par
# (choose mpi_size to divide nFilterCycles)
srun -n <mpi_size> -c <nThreads> --cpu-bind=cores ./Filter_mpi.x > run.dat 2> error.dat
```

Outputs quasiparticle energies (`eval.dat`), eigenvectors (`psi.dat`), and a
monolithic `output.dat` consumed by the BSE step.

### 2. BSE (`bse_cplx`)

```sh
source LOAD_COMPILE_ENV_NVHPC.sh    # PrgEnv-nvidia + cray-fftw + cray-libsci
cd bse_cplx && make                 # produces bse_cplx_gpu.x (GPU + ScaLAPACK)
# or: make -f Makefile_mkl          # produces bse_cplx.x (CPU-only, single node)

# Run in a directory holding the filter output + input.par (use >= 2 ranks)
export BSE_GPU_MEM_GB=36
srun -n 4 --gpus-per-task=1 --gpu-bind=closest -c 32 ./bse_cplx_gpu.x > run.dat 2> error.dat
```

Outputs the exciton spectrum (`exciton.dat`), oscillator strengths (`OS.dat`,
`M.dat`, `rs.dat`), and angular-momentum character.

See the [documentation](https://dchabeda.github.io/FilterBSE/) for the full input
reference, job options, parallelization guide, and output formats. Worked examples
live under `filter_mpi/EXAMPLES/`.

## Requirements

- C compiler with MPI and OpenMP (Cray `cc` wrapper on Perlmutter, over Intel or
  NVHPC)
- FFTW 3
- LAPACK/BLAS: Intel MKL, or Cray libsci (for the ScaLAPACK / GPU `bse_cplx` build)
- NVHPC + an NVIDIA GPU (cc80 / A100) for the optional `bse_cplx` GPU offload

## Documentation

The docs are built with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/)
from the `docs/` folder and published automatically to GitHub Pages on every push
to `main`. To build them locally:

```sh
pip install -r docs/requirements.txt
mkdocs serve      # live preview at http://127.0.0.1:8000
```

## License

Released under the [MIT License](LICENSE). Copyright © 2026 Daniel Chabeda.
