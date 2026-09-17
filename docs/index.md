# FilterBSE

FilterBSE is a suite of electronic-structure codes for semiconductor
nanostructures. It computes single-particle states by **filter
diagonalization** and optical / excitonic properties by solving the
**Bethe–Salpeter equation (BSE)**.

The pipeline runs in two steps:

1. **Filter** — compute quasiparticle states and energies on a real-space grid.
2. **BSE** — build the electron–hole kernel from those states and diagonalize it
   for excitons and optical spectra.

## Where to start

<div class="grid cards" markdown>

- :material-cube-outline: **[filter_mpi](filter_mpi/index.md)**

    The MPI filter-diagonalization code. Distributed filter → orthogonalize →
    diagonalize pipeline that produces quasiparticle states (`psi.dat`) and
    energies (`eval.dat`).

- :material-atom-variant: **bse_cplx** *(coming next)*

    The complex Bethe–Salpeter solver that consumes the filter output to compute
    excitons and optical properties.

</div>

## Repository map

| Directory | What it is |
|---|---|
| `filter_mpi/` | MPI filter-diagonalization (documented here first) |
| `filter/` | Serial / reference filter code |
| `periodic_filter/` | Periodic-system filter |
| `bse_cplx/` | Complex BSE solver (docs next) |
| `bse_real/`, `bse_gpu/` | Real-valued and GPU BSE variants |
| `pots/` | Pseudopotential data |
| `auxiliary/` | Helper tools |

## Building and running (quick orientation)

The codes are C + MPI, built with per-directory Makefiles. On Perlmutter the
build environment is loaded with:

```sh
source LOAD_COMPILE_ENV.sh   # module load PrgEnv-intel, cray-fftw
```

See **[filter_mpi → Job options](filter_mpi/job-options.md)** for how to build,
launch, and configure a run, and **[Input files](filter_mpi/input-files.md)** for
the input format.
