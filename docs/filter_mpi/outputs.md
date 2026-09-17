# Outputs

This page lists the files `filter_mpi` writes. The three that matter most are
**`eval.dat`** (energies), **`psi.dat`** (eigenvectors), and **`output.dat`**
(the monolithic file the BSE codes read).

## Primary results

### `eval.dat` — eigenvalues

Text, one line per state:

```
index    eigenvalue(%.16lg)    sigma_E(%lg)
```

`sigma_E` is the eigenvalue variance (ghost-state diagnostic). Periodic runs
write one file per k-point: `eval-<k>.dat`.

### `psi.dat` — eigenvectors

Raw binary: a contiguous `fwrite` of the grid-basis eigenvectors,
`mn_states × complex_idx × nspinngrid` doubles, one state after another. Periodic
runs write `psi-<k>.dat`.

!!! info "Byte layout is a contract"
    State \(g\) begins at byte offset \(g\cdot\texttt{stlen}\cdot 8\) with
    \(\texttt{stlen} = \texttt{complex\_idx}\cdot\texttt{nspinngrid}\). The
    serial writer, the distributed MPI-IO writer (`dist_write_psi_dat`), the
    `get_n_states.x` utility, and the BSE readers all agree on this layout. If
    you post-process `psi.dat` by hand, honor it exactly.

### `output.dat` — monolithic BSE input

Written by `save_output` when `saveOutput != 0`. It packs everything the BSE
codes need into one file: an `output_tag`, the `ist` sizes (`ngrid`,
`nspinngrid`, `mn_states_tot`, `natoms`, `n_atom_types`, atom types, `nspin`,
`complex_idx`), selected `par` values (`KE_max`, `fermi_E`), flags (SO, NL, LR,
`useSpinors`, `isComplex`), atomic coordinates, full grid parameters and the
`fwrite` of the grid x/y/z, then `fwrite` of `eig_vals`, `sigma_E`, and the
eigenvectors, terminated by an `EOF` marker.

!!! warning "Not produced on the distributed (MPIOrtho) path"
    The `MPIOrtho = 1` path does **not** write `output.dat`. Instead it prints a
    note to run the BSE code from `psi.dat` / `eval.dat` / `conf.par` (the
    BSE "unsafe" input path). Plan for this when running at scale.

## Restart / intermediate binaries

| File | Contents |
|---|---|
| `psi-filt-<jns>-<rank>.dat` | per-rank filtered states (stage 1) |
| `psi-filt.dat` / `psi-filt-<k>.dat` | filtered states for restart |
| `psi-diag.dat` | diagonalized states (for `restartFromSigma`) |
| `checkpoint_1/2/3.dat` | full job state between stages (`saveCheckpoints = 1`) |

See [Job options → Restarting](job-options.md#restarting-a-run).

## Cube files

Gaussian cube format, written when the relevant flags are set:

- `local-pot.cube` — the local pseudopotential on the grid (from `mod_pseudopot`).
- Per-state wavefunction cubes via `write_cube_file` / `print_bloch_cubes_k` when
  `printCubes = 1` (the second token sets how many). The reference example
  produces e.g. `homo-0-Up.cube`, `lumo+0-Up.cube`.

## Diagnostics

Various `.dat` files written along the way (useful for debugging and
convergence checks):

- `ene-filt-jns-<jns>-<rank>.dat` — filtered-state energies vs targets
- `ene_targets.dat` — the energy targets
- `zn.dat` — Chebyshev support points; `func.dat` — the filter function
- `Emin-init.dat`, `Emax-init.dat`, `eval_aux.dat` — energy-range diagnostics
- `grid.dat`, `ksqr.dat` — grid and \(k^2\)
- `G_vecs.dat`, `list_NL_grid.dat` — periodic / non-local grids
- `conf.dat` — echoed atomic configuration
- `hmat.dat`, `projectors.dat`, `strain.dat`, `angular.dat` — stage-specific dumps

## stdout / stderr

Rank 0 prints a banner and a numbered stage log:

```
1. INITIALIZING JOB
2. CALCULATING HAMILTONIAN ENERGY RANGE
3. GENERATING COEFFICIENTS
4. RUN FILTER CYCLE        (with progress bars over Chebyshev iterations)
5. ORTHOGONALIZATING
6. DIAGONALIZING HAMILTONIAN
7. CALCULATING VARIANCE OF EIGENVALUES
POST-PROCESSING
DONE                       (with CPU and wall times)
```

At startup `print_input_state` echoes the parsed job configuration so the log is
self-documenting. Errors go to `stderr` and typically trigger `MPI_Abort` /
`exit(EXIT_FAILURE)`. In the example runs stdout is redirected to `run.dat` and
stderr to `error.dat`.
