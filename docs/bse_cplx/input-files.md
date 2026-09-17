# Input files

`bse_cplx` always reads a parameter file **`input.par`**. The quasiparticle
wavefunctions and metadata come from the filter output, through one of two paths
selected by the first command-line argument.

## Two input paths

The first CLI argument sets `initUnsafe`:

```sh
./bse_cplx_gpu.x        # initUnsafe = 0  → "safe" path (default)
./bse_cplx_gpu.x 1      # initUnsafe = 1  → "unsafe" path
```

=== "Safe path (default)"

    Reads the **monolithic `output.dat`** written by `filter_mpi` (`saveOutput`).
    The header — grid, atoms, `eig_vals`, `sigma_E`, and the flags/sizes — is read
    up front; the code records the byte offset of the wavefunction block and
    verifies a trailing `EOF` marker. The wavefunctions themselves are **not**
    read yet — only the selected QP states are loaded later, on demand, by
    `get_qp_basis` seeking into `output.dat`.

=== "Unsafe path (`argv[1] = 1`)"

    Reads separate raw files instead of `output.dat`:

    | File | Contents |
    |---|---|
    | `unsafe_input.par` | grid + sizes keyword file (see below) |
    | `eval.par` | ASCII `idx eigval sigma` per state |
    | `conf.par` | atomic configuration `symbol x y z` |
    | `psi.par` | raw binary wavefunctions (`double complex`) |

    This is the path to use with the [`filter_mpi` distributed
    (MPIOrtho)](../filter_mpi/parallelization.md) output, which produces
    `psi.par` / `eval.par` / `conf.par` rather than a monolithic `output.dat`.

    !!! warning "Use at your own risk"
        The code itself flags that in the unsafe path the wavefunctions "might
        not be aligned with the grid." `psi.par` must be exactly
        `mn_states_tot × nspinngrid × sizeof(double complex)` bytes or the run
        aborts. Note the `nHoles` / `nElecs` values in `unsafe_input.par` are
        **hints only** — the real hole/electron window is chosen later by the
        \(\sigma_E\)/Fermi logic, which correctly skips interior ghost states a
        contiguous window could not.

## `input.par` reference

Every line is `key = value`. An unknown key prints the full list of allowed keys
and exits. A missing `input.par` aborts.

Here is the example from `EXAMPLES/1x1x1_CsPbI3/bse`:

```ini
maxHoleStates = 10
maxElecStates = 10
sigmaECut = 0.0001
epsX = 6.1
epsY = 6.1
epsZ = 6.1
nThreads = 1
spinOrbit = 1
NonLocal = 1
calcDarkStates = 0
calcSpinAngStat = 0
timingSpecs = 1
fermiEnergy = -0.19
deltaEhole = 0.04
deltaEelec = 0.04
KEmax = 10.0
coulombDone = 0
restartCoulomb = 0
gpuAccel = 1
```

### QP basis window

| Key | Type | Meaning | Default |
|---|---|---|---|
| `fermiEnergy` | double | energy dividing holes from electrons (a.u.) | -0.18 |
| `sigmaECut` | double | max eigenvalue variance \(\sigma_E\) to accept a state (rejects ghosts) | 0.0001 |
| `maxHoleStates` | long | cap on hole states nearest the gap (−1 = unconstrained) | -1 |
| `maxElecStates` | long | cap on electron states nearest the gap (−1 = unconstrained) | -1 |
| `deltaEhole` | double | desired hole-band energy span (advisory only) | 0 |
| `deltaEelec` | double | desired electron-band energy span (advisory only) | 0 |

!!! note "How the window works"
    A state is a **hole** if `sigma_E < sigmaECut && eigval < fermiEnergy`, an
    **electron** if `sigma_E < sigmaECut && eigval > fermiEnergy`. The \(\sigma_E\)
    cut simultaneously discards unconverged/ghost filter states. The exciton
    basis size is `n_xton = n_holes × n_elecs`.

### Screening & interaction

| Key | Type | Meaning | Default |
|---|---|---|---|
| `epsX`, `epsY`, `epsZ` | double | anisotropic dielectric constants for the screened (direct) Coulomb term | 0 (must be ≥ 0) |
| `KEmax` | double | maximum kinetic energy | 20.0 |
| `longRange` | int | potentials include long-range terms (no truncation) | 0 |
| `spinOrbit` | int | spin–orbit (forces `useSpinors=1`, `NonLocal=1`) | 0 |
| `useSpinors` | int | 2-component spinor wavefunctions (→ complex, `nspin=2`) | 0 |
| `NonLocal` | int | non-local potential | 0 |

### Output / analysis control

| Key | Type | Meaning | Default |
|---|---|---|---|
| `calcSpinAngStat` | int | compute spin / angular-momentum statistics | 1 |
| `calcDarkStates` | int | compute the dark (spin-forbidden) manifold | 0 |
| `printFPDensity` | int | print fixed-point exciton densities | 0 |
| `timingSpecs` | int | print per-step BSE transform timings | 0 |

### Parallelization

| Key | Type | Meaning | Default |
|---|---|---|---|
| `nThreads` | long | OpenMP threads per rank (FFTW, dipole/angular loops, threaded ScaLAPACK) | — |
| `gpuAccel` | int | offload the Coulomb kernel to GPU when available (0 forces CPU) | 1 |

See [Parallelization](parallelization.md) for how these combine with rank counts
and `BSE_GPU_MEM_GB`.

### Run mode / restart

| Key | Type | Meaning | Default |
|---|---|---|---|
| `noCalcExciton` | int | stop after single-particle (QP) properties | 0 |
| `calcCoulombOnly` | int | exit after computing the kernel | 0 |
| `coulombDone` | int | kernel already computed → load from disk, skip to BSE | 0 |
| `restartCoulomb` | int | resume the kernel from per-rank dumps | 0 |
| `saveCheckpoints` | int | save job checkpoints | 0 |
| `restartFromChk` | int | checkpoint ID to restart from | -1 |

### `unsafe_input.par` keywords

Used only on the unsafe path: `mnStatesTot`, `nHoles`, `nElecs`, `nAtoms`, `nx`,
`ny`, `nz`, `xmin`/`ymin`/`zmin`, `dx`/`dy`/`dz`, `fermiE`, `sigmaECut`.
(`nHoles`/`nElecs` are hints only.)

## Kernel restart files

Produced by earlier runs and re-read to skip recomputing the Coulomb kernel (see
[Outputs](outputs.md) and [Job options → Restarting](job-options.md#restarting-and-checkpointing)):

- `direct.dat` / `exchange.dat` — whole-matrix reload (`coulombDone = 1`)
- `direct-<rank>.dat` / `exchange-<rank>.dat` — per-rank dumps (`restartCoulomb = 1`)
