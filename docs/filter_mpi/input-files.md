# Input files

A `filter_mpi` run reads its configuration from a small set of plain-text files
in the working directory. The main one is **`input.par`**; the atomic geometry is
in **`conf.par`**; the pseudopotentials come from **`pot<Sym>.par`** and (for
spin–orbit / non-local) **`SO_<Sym>.par`** files. Periodic runs additionally read
**`periodic_input.par`**.

!!! warning "`input.par` must exist"
    The program aborts immediately if `input.par` is not found in the working
    directory. `conf.par` is likewise required.

## `input.par` format

Every line is `key = value` (case-sensitive), parsed with `fscanf(pf, "%s = %s", ...)`.
A handful of keys read **extra tokens on the same line** after the value (noted
below). An **unknown key prints the full list of allowed keys and exits**, so a
typo is caught immediately.

Here is a complete, working non-periodic example (from
`EXAMPLES/1x1x1_CsPbI3`):

```ini
nx = 44
ny = 44
nz = 44
dGrid = 0.5
mStatesPerFilter = 16
nFilterCycles = 16
nCheby = 2048
VBmin = -0.50
VBmax = -0.20
CBmin = -0.167
CBmax = -0.10
nThreads = 16
hamThreads = 1
spinOrbit = 1
NonLocal = 1
interpolatePot = 0
setTargets = 1 10 6
calcPotOverlap = 0
getAllStates = 1
sigmaECut = 0.001
timeHamiltonian = 1
fermiEnergy = -0.19
setSeed = 0
printNorm = 0
printCubes = 1 3
saveCheckpoints = 0
```

!!! tip "Energies are in atomic units (Hartree)"
    `VBmin/VBmax/CBmin/CBmax`, `fermiEnergy`, and `KEmax` are all in a.u. Grid
    spacings (`dGrid`, `dx`…) are in Bohr.

### Grid and geometry

| Key | Type | Meaning | Default |
|---|---|---|---|
| `nx`, `ny`, `nz` | long | grid points per axis (a starting estimate — the code shrinks the grid to fit the nanocrystal) | required |
| `dGrid` | double | uniform grid spacing (sets `dx=dy=dz`, Bohr) | — |
| `dx`, `dy`, `dz` | double | per-axis spacing (Bohr) | from `dGrid` |
| `box_z` / `box_Z` | double | periodic box dimension (Bohr) | 0.0 |
| `centerConf` | int | center the structure at the center of mass | 1 |
| `periodic` | int | 0 = 0D cluster, 1 = periodic | 0 |
| `diagKIdx` | int | k-point index to diagonalize | 0 |

### Filter algorithm

| Key | Type | Meaning | Default |
|---|---|---|---|
| `nFilterCycles` | long | number of random initial states (filter cycles) | required |
| `mStatesPerFilter` | long | number of energy targets per filter | required |
| `nCheby` | long | number of Chebyshev/Newton terms | required |
| `VBmin`, `VBmax` | double | valence-band energy window (a.u.) | — |
| `CBmin`, `CBmax` | double | conduction-band energy window (a.u.) | — |
| `KEmax` | double | maximum kinetic energy | 20.0 |
| `fermiEnergy` | double | Fermi energy (a.u.) | -0.18 |
| `setTargets` | int + `nVB nCB` | manually set the VB/CB target counts (must sum to `mStatesPerFilter`) | 0 |
| `approxEnergyRange` | int | approximate the spectral bound from the local potential | 0 |
| `setSeed` | int + `seed` | fix the RNG seed for reproducibility | 0 |
| `inputPsiFilt` | int + `start end` | read initial states from disk instead of random | 0 |

!!! note "Placing targets"
    It wastes effort to place filters inside the gap: set `VBmax` just above the
    HOMO and `CBmin` just below the LUMO. Because the valence band has a higher
    density of states, put more targets in the VB — e.g. `setTargets = 1 10 6`
    (10 VB + 6 CB = 16 = `mStatesPerFilter`).

### Pseudopotential

| Key | Type | Meaning | Default |
|---|---|---|---|
| `useStrain` | int | strain-dependent pseudopotentials | 0 |
| `interpolatePot` | int | interpolate cubic/ortho potentials (uses `_cubic`/`_ortho` files) | 0 |
| `longRange` | int | potentials include a long-range term (no truncation) | 0 |
| `localStructureDependent` | int | local-structure-dependent (Δv) correction | 0 |
| `scaleSurfaceCs` | double | fractional scaling of surface Cs | 1.0 |
| `crystalStructure` | string | e.g. `wurtzite`, `zincblende` | unknown |
| `outmostMaterial` | string | e.g. `CdSe`, `InP` | unknown |
| `readProj` | int | read projector data from files | 0 |
| `readCoeffs` | int | read Newton coefficients / sample points from files | 0 |
| `psiZeroCut` | double | wavefunction zero cutoff | 1e-16 |

### Spin / interaction

| Key | Type | Meaning | Default |
|---|---|---|---|
| `useSpinors` | int | use 2-component spinors (→ `nspin=2`, complex) | 0 |
| `spinOrbit` | int | spin–orbit coupling (forces `useSpinors=1`) | 0 |
| `NonLocal` | int | non-local (Kleinman–Bylander) potential | 0 |
| `noTimeRev` | int | 1 = do **not** apply time reversal | 1 |

!!! info "Derived flags"
    `spinOrbit=1` forces `useSpinors=1`; `useSpinors=1` (or `periodic=1`) forces
    a complex calculation (`isComplex=1`, `complex_idx=2`). With time reversal
    active (`useSpinors && !noTimeRev`) the basis is doubled.

### Parallelization

| Key | Type | Meaning | Default |
|---|---|---|---|
| `nThreads` | long | total OpenMP threads per rank | — |
| `hamThreads` | int | OpenMP/FFTW threads for the Hamiltonian | 1 |
| `nestedOMP` | int | enable nested OpenMP | 0 |
| `useFastHam` | int | fast Hamiltonian path | 0 |
| `useMPIOMP` | int | hybrid MPI+OMP mode | 0 |
| `MPIOrtho` | int | 0 = single-node ortho, 1 = distributed ortho | 0 |
| `MPIDiag` | int | 0 = single-node diag, 1 = distributed diag | 0 |
| `parallelSigma` | int | 1 = parallel over states, 0 = parallel Hamiltonian | 1 |

See [Parallelization](parallelization.md) for how these interact and how to
choose rank counts.

### Output control

| Key | Type | Meaning | Default |
|---|---|---|---|
| `printPsiFilt` | int | write filtered states | 1 |
| `printPsiOrtho` | int | write orthogonalized states | 0 |
| `printPsiDiag` | int | write diagonalized states | 0 |
| `printCubes` | int + `ncubes` | write cube files | — |
| `printGaussCubes` | int + `start end` | write Gaussian-basis cubes | — |
| `printNorm` | int | print norms every 100 Chebyshev iterations | 0 |
| `calcSpinAngStat` | int | spin / angular-momentum statistics (SO only) | 0 |
| `saveOutput` | int | write the monolithic `output.dat` | 1 |
| `getAllStates` | int + optional `sigma_E_cut` | write all states vs only converged ones | 1 |
| `sigmaECut` | double | variance cutoff defining a "converged" state | 0.01 |
| `calcPotOverlap` | int | compute \(\langle i|V|j\rangle\) | 0 |
| `timeHamiltonian` | int | print Hamiltonian timing | 0 |
| `fftWisdomDir` | string | directory for `fftw_wisdom.dat` | "" |

### Restart / checkpoint

| Key | Type | Meaning | Default |
|---|---|---|---|
| `saveCheckpoints` | int | write `checkpoint_1/2/3.dat` between stages | 0 |
| `restartFromCheckpoint` | int | checkpoint ID (0–4) to resume from | 0 |
| `restartFromOrtho` | int + `nStates` | read `psi-filt.dat`, start at ortho | 0 |
| `restartFromSigma` | int + `nStates` | read `psi-diag.dat`, start at sigma | 0 |
| `retryFilter` | int | retry (adds 16 cycles, 1024 cheby) if 0 states found | 0 |

See [Job options → Restarting](job-options.md#restarting-a-run) for how these map
onto the pipeline stages.

## `conf.par` — atomic configuration

First line is the atom count; each following line is `symbol x y z`:

```
15
Cs 0.000000 11.884000 -5.942000
I  5.942000 11.884000  0.000000
I  5.942000  5.942000 -5.942000
...
Pb 5.942000  5.942000  0.000000
```

Coordinates are Cartesian. This file is required.

## Pseudopotential files

The local radial pseudopotential for each element `<Sym>` is read from
`pot<Sym>.par`, a two-column `r  value` table. Depending on `input.par` flags the
code looks for variants:

| File | When used |
|---|---|
| `pot<Sym>.par` | default local potential |
| `pot<Sym>_SR.par`, `pot<Sym>_LR.par` | short-/long-range split (`longRange`) |
| `pot<Sym>_a4.par`, `pot<Sym>_a5.par` | strain parameters (`useStrain`) |
| `pot<Sym>_cubic.par`, `pot<Sym>_ortho.par` (+ `_SR`/`_LR`) | interpolation (`interpolatePot`) |

## Spin–orbit / non-local files

When `spinOrbit` or `NonLocal` is on, each non-ligand element reads
`SO_<Sym>.par`, containing three values: the SO parameter and two non-local
parameters. Interpolation variants `SO_<Sym>_cubic.par` etc. are used with
`interpolatePot`. Ligand atoms (e.g. passivants) get no SO/NL terms.

The `EXAMPLES/1x1x1_CsPbI3/filter` directory shows a full set: `potCs.par`,
`potI.par`, `potPb.par`, `SO_Cs.par`, `SO_I.par`, `SO_Pb.par` plus `_cubic`
variants.

## Periodic inputs (`periodic = 1`)

Read only when `periodic = 1`:

- **`periodic_input.par`** — lattice and k-mesh. Keys: `a`, `b`, `c`, `alpha`,
  `beta`, `gamma`; lattice vectors `a1`, `a2`, `a3` and reciprocal vectors `b1`,
  `b2`, `b3`; `nBands`, `nbMin`, `nbMax`; the k-mesh `nk1`, `nk2`, `nk3`;
  `readKPath`.
- **`kpath.par`** — an explicit k-path, read when `readKPath = 1`.

## Restart binaries

These are produced by earlier runs and re-read on restart (see
[Outputs](outputs.md) and [Job options](job-options.md#restarting-a-run)):

- `psi-filt.dat` / `psi-filt-<k>.dat` — filtered states
- `psi-diag.dat` — diagonalized states
- `checkpoint_1/2/3.dat` — full job state between stages
- `zn.dat` + coefficient files — read when `readCoeffs = 1`
