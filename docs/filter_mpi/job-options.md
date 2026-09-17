# Job options

This page covers building the code, launching it, and the run modes you select
through `input.par`. All configuration is via the input files — **the production
binary reads no command-line flags and no environment variables** (only MPI's own
and the OpenMP `OMP_*` variables affect it).

## Building

The build environment on Perlmutter is Intel + Cray FFTW:

```sh
source LOAD_COMPILE_ENV.sh
# module load PrgEnv-intel/8.5.0
# module load cray-fftw/3.3.10.8
```

Then, in `filter_mpi/`:

```sh
make            # builds Filter_mpi.x (default target)
make cube       # builds makecube.x (post-processing cube generator)
make clean      # removes *.o *.x *.d
```

Build details (from `Makefile`):

- Compiler `cc` (the Cray wrapper over Intel).
- Flags: `-DMKL_ILP64 -DFFTW_FFT -O3 -qopenmp -mavx2` (plus `-MMD -MP` for
  auto-generated header dependencies).
- Links MKL (ILP64), OpenMP, pthreads, and threaded FFTW
  (`-lfftw3_omp -lfftw3`).

!!! note "get_n_states.x is separate"
    The utility `get_n_states.x` has its own `main` and is **not** built by the
    Makefile — compile it separately. It slices states out of a `psi.dat`-style
    binary: `get_n_states start end nspinngrid cmplx filename`.

## Launching

`filter_mpi` uses **state (column) parallelism**: each rank filters
`nFilterCycles / mpi_size` of the random starting vectors, holding whole
grid-length state vectors (the grid is *not* split across ranks). Within a rank,
OpenMP (`nThreads`) and threaded FFTW (`hamThreads`) accelerate the Hamiltonian.

A representative launch (adapt the SLURM header to your allocation):

```sh
export OMP_NUM_THREADS=16
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

srun -n <mpi_size> -c <nThreads> --cpu-bind=cores \
     ./Filter_mpi.x > run.dat 2> error.dat
```

!!! warning "Pick `mpi_size` to divide `nFilterCycles`"
    `n_filters_per_rank = nFilterCycles / mpi_size`. If `mpi_size` does not
    divide `nFilterCycles` (or exceeds it), ranks fall back to 1 filter each and
    you waste resources. Also set `nThreads × mpi_size ≈ cores per node`. For
    **periodic** jobs, `mpi_size = n_k_pts × (a divisor of nFilterCycles)`, and
    `mpi_size` must be a multiple of `n_k_pts`. See
    [Parallelization](parallelization.md) for the full story.

A typical run directory contains: `input.par`, `conf.par`, the `pot<Sym>.par` and
`SO_<Sym>.par` files, and (if periodic) `periodic_input.par`.

## Run modes

The pipeline is driven by a fall-through switch on `restartFromCheckpoint`
(0 → 4) in `main.c`. The flags below select which path runs.

### Cluster vs periodic

`periodic` selects the entire k-aware pipeline: G-vectors, a k-point mesh,
k-communicators (`MPI_Comm_split`), a Bloch Hamiltonian, and per-k output
(`eval-<k>.dat`, `psi-<k>.dat`). With `periodic = 0` the code runs the 0D cluster
path.

### Orthogonalization / diagonalization path

| Flags | Path taken |
|---|---|
| default (`MPIOrtho=0`, `MPIDiag=0`) | Gather all filtered states to rank 0, then serial SVD ortho (`zgesvd`/`dgesvd`) and serial diag (`zheev`/`dsyev`). |
| `MPIOrtho = 1` | Fully distributed ortho **and** diag (`dist_linalg.c`): states are never gathered; ring/Gram SVD, distributed diagonalization, MPI-IO output. |
| `MPIDiag = 1` | Distributed subspace-Hamiltonian construction (`diag_H_mpi`), diagonalized on one rank. |

!!! tip "When to use `MPIOrtho = 1`"
    Use the distributed path when the full set of filtered states would not fit
    in one node's RAM (the serial gather path caps out around 450 GiB/node). Be
    aware it uses a coarser SVD cutoff (`1e-7` vs the serial `1e-10`), so it
    retains slightly fewer states — see the
    [Parallelization page](parallelization.md#stage-2-distributed-orthogonalization-svd-via-the-gram-matrix).
    On an MPIOrtho restart, `nStates` must be divisible by `mpi_size`.

### Filter-only

Set `calcFilterOnly = 1` to stop after the filter stage (useful for staging large
runs or debugging the filter).

### Eigenvalue variance (sigma)

After diagonalization the code computes each eigenvalue's variance
\(\sigma_E = \sqrt{\langle\psi|\hat H^2|\psi\rangle - \langle\psi|\hat H|\psi\rangle^2}\)
as a ghost-state / convergence diagnostic. `parallelSigma` chooses between
parallel-over-states (default) and parallel-Hamiltonian evaluation depending on
memory. States with \(\sigma_E \le\) `sigmaECut` are considered converged.

## Restarting a run

Restarts let you resume at a pipeline boundary instead of re-filtering:

| Setting | Reads | Resumes at |
|---|---|---|
| `restartFromOrtho = 1 <nStates>` | `psi-filt.dat` | orthogonalization (forces checkpoint 1) |
| `restartFromSigma = 1 <nStates>` | `psi-diag.dat` | the sigma stage (forces checkpoint 3) |
| `restartFromCheckpoint = 4` | eval/psi from disk | post-processing / output only |
| `saveCheckpoints = 1` | — | writes `checkpoint_1/2/3.dat` between stages so any of the above can resume |

!!! note
    Checkpointing is not supported on the periodic k-grouped path except via
    `restartFromOrtho`. `retryFilter = 1` will automatically add 16 filter cycles
    and 1024 Chebyshev terms and rerun if the filter finds zero eigenstates.
