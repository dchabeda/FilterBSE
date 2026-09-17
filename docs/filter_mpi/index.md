# filter_mpi

`filter_mpi` is the MPI + OpenMP filter-diagonalization code. It computes
single-particle (quasiparticle) states and energies of a nanostructure on a
real-space grid using a semiempirical pseudopotential Hamiltonian, and writes
output that the downstream BSE codes consume.

It runs a three-stage pipeline:

<div class="grid cards" markdown>

- :material-filter: **1. Filter**

    A Chebyshev/Newton filter projects random starting vectors onto a set of
    energy targets in the valence/conduction windows, producing a redundant set
    of filtered states.

- :material-vector-arrange-below: **2. Orthogonalize**

    An SVD (serial) or distributed Gram-matrix SVD removes the numerical
    redundancy, leaving an orthonormal basis.

- :material-grid: **3. Diagonalize**

    The Hamiltonian is built and diagonalized in that small basis, giving
    quasiparticle energies and grid-basis eigenvectors, plus an eigenvalue
    variance \(\sigma_E\) ghost-state check.

</div>

## How to read these docs

| Page | What's in it |
|---|---|
| [Input files](input-files.md) | `input.par` keyword reference, `conf.par`, pseudopotential/SO files, periodic inputs |
| [Job options](job-options.md) | Building, launching, run modes, restart, path selection |
| [Codebase & module layout](codebase.md) | Every source file, `main.c` control flow, headers |
| [Pipeline: filter → ortho → diag](pipeline.md) | What each stage does and how data flows between them |
| [Parallelization](parallelization.md) | MPI decomposition, the distributed ring/Gram scheme, rank-count tips |
| [Outputs](outputs.md) | Every file the code writes and its format |

## Quick start

```sh
# 1. Load the build environment (Perlmutter)
source LOAD_COMPILE_ENV.sh          # PrgEnv-intel + cray-fftw

# 2. Build
cd filter_mpi && make               # produces Filter_mpi.x

# 3. Prepare a run directory with input.par, conf.par,
#    pot<Sym>.par and SO_<Sym>.par files (see Input files)

# 4. Launch with MPI (mpi_size should divide nFilterCycles)
srun -n <mpi_size> -c <nThreads> --cpu-bind=cores ./Filter_mpi.x > run.dat 2> error.dat
```

See [Job options](job-options.md) for the details behind each step.
