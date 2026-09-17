# Codebase & module layout

`filter_mpi` is C + MPI + OpenMP, using MKL (LAPACK/BLAS, ILP64) and FFTW. The
build produces one production binary, **`Filter_mpi.x`**, plus two standalone
utilities.

## Top-level control flow (`main.c`)

`main()` is the sole production entry point. Its flow:

1. `MPI_Init`, get rank/size, initialize default k-group fields.
2. **`mod_init`** — read input, atoms, grid, periodic setup, energy targets.
3. **`mod_mem_alloc`** — allocate the psi arrays, coefficients, eigenvalues.
4. **`mod_pseudopot`** — build the local potential and NL/SO projectors.
5. A **`switch (restartFromCheckpoint)`** with fall-through cases 0 → 4 that
   drives the pipeline:

    | Case | Stage | What happens |
    |---|---|---|
    | 0 | filter | `mod_filter`, then branch: `calcFilterOnly` exits; `periodic` → `run_periodic_postfilter`; `MPIOrtho` → `run_dist_postfilter`; else gather to rank 0 (+ optional checkpoint). |
    | 1 | ortho | serial `restart_from_ortho` + `mod_ortho` (periodic / MPIOrtho variants for those paths). |
    | 2 | diag | `mod_diag`. |
    | 3 | sigma | rank-0 `mod_sigma` (skipped when periodic — variance is done inside `mod_diag`). |
    | 4 | output | rank-0 `mod_output`, or load-and-output post-processing. |

6. Optional `mod_optional_output` on rank 0 (cube files, diagnostics).
7. Free memory, print timing, `MPI_Finalize`.

The cases **fall through**, so a fresh run (case 0) proceeds through ortho → diag
→ sigma → output in one invocation; restarts jump in partway.

## Module drivers (`mod_*`)

These are the pipeline stages, orchestrated by `main.c`:

| File | Role |
|---|---|
| `mod_init.c` | read input → read conf → build grid → periodic init → k-communicators → energy targets |
| `mod_mem.c` | allocate psi/phi/potential/projector/eigenvalue arrays |
| `mod_pseudopot.c` | build local pseudopotential on the grid + projectors; writes `local-pot.cube` |
| `mod_filter.c` | **Stage 1** — energy range → Chebyshev/Newton coeffs → random init states → filter cycles |
| `mod_ortho.c` | **Stage 2** — single-rank SVD orthogonalization |
| `mod_portho.c` | block/distributed ortho (largely superseded by `dist_linalg.c`) + `read_psi_from_disk` |
| `mod_diag.c` | **Stage 3** — subspace diagonalization; also the periodic post-filter pipelines |
| `mod_sigma.c` | eigenvalue variance (ghost-state check); state-parallel vs H-parallel |
| `mod_output.c` | writes `eval.dat` / `psi.dat` / `output.dat`; `select_conv_states` |
| `mod_optional_output.c` | cube files, extra diagnostics, post-processing output |
| `dist_linalg.c` | the **distributed (MPIOrtho) pipeline**: ring Gram matrix, distributed SVD ortho, distributed diag, MPI-IO `psi.dat`, `run_dist_postfilter` |

## Compute & support files

| File | Role |
|---|---|
| `read.c` | all input parsing (`read_input`, `read_conf`, `read_pot`, `read_periodic_input`) |
| `init.c` | grid params, `ksqr`, energy targets, random `init_psi`; writes `grid.dat`, `ksqr.dat` |
| `init_periodic.c` | G-vectors, k-points, `setup_k_communicators` |
| `filter.c` | Chebyshev filter cycles (`run_filter_cycles`, `filter_cycle`, `p_hnorm`); MPI gathers |
| `energy.c` | energy of filtered states; `get_energy_range` (spectral bounds) |
| `coeff.c` | Newton interpolation coefficients for the filter function; writes `zn.dat` |
| `hamiltonian.c` | non-periodic \(\hat H\) application (kinetic via FFT + local + NL + SO) |
| `phamiltonian.c` | periodic (Bloch, k-dependent) Hamiltonian |
| `ghamiltonian.c` | Gaussian-basis / auxiliary Hamiltonian |
| `ortho.c` | serial SVD ortho (`ortho_real`/`ortho_cplx`), cutoff at `SVDEPS = 1e-10` |
| `Hmat.c` | serial subspace-Hamiltonian build + `diag_H` (`dsyev`/`zheev`) |
| `Hmat_mpi.c` | `diag_H_mpi` — distributed subspace-Hamiltonian construction |
| `norm.c` | `calc_norm`, `normalize`, `normalize_all` |
| `save.c` | checkpoints, `print_input_state`, `save_output` (writes `output.dat`) |
| `write.c` | file writers (`write_eval_dat`, `write_psi_dat`, `write_cube_file`, …) |
| `projectors.c` | SO/NL projector generation |
| `strain.c` | strain-dependent pseudopotential corrections |
| `interpolate.c` | radial potential interpolation |
| `rand.c` | RNG (`ran_nrc`, `Randomize`) |
| `aux.c` | timing, progress bars, small helper diagonalization |
| `vector.c` | 3-vector math |
| `integral.c` | Gaussian-basis overlap/kinetic/potential integrals |
| `nerror.c` | `nerror()` fatal-error helper |

## Headers

- **`fd.h`** — the master header (`#pragma once`): all library includes, **every
  struct definition** (`flag_st`, `index_st`/`ist`, `par_st`, `grid_st`,
  `atom_info`, `pot_st`, `zomplex`, `parallel_st`, …), macros
  (`AUTOEV`, `ANGTOBOHR`, `SVDEPS 1e-10`, `N_MAX_ATOM_TYPES 20`, complex
  arithmetic), and most function prototypes.
- **`main.h`** — includes `fd.h` plus every module header.

## Standalone utilities (not in `Filter_mpi.x`)

| File | Binary | Purpose |
|---|---|---|
| `makecube.c` (+ `size.c`) | `makecube.x` (`make cube`) | post-processing cube generator |
| `get_n_states.c` | `get_n_states.x` (built separately) | slice states out of a `psi.dat` binary: `get_n_states start end nspinngrid cmplx filename` |
| `dist_linalg.c` (with `-DDIST_LINALG_TEST`) | self-test `main` | unit test for the ring Gram/SVD (run on N ranks) |

## Legacy / experimental (present but not linked)

The Gaussian-basis path and some BSE-projection helpers are present but **not**
compiled into `Filter_mpi.x` (the `gauss_driver` call in `main.c` is commented
out). Treat these as experimental unless verified:

`gauss.c`, `dipole.c`, `pot_mat_elems.c`, `projection_mat.c`,
`pseudopotential.c`, `angularTest.c`. Also, the post-ortho re-diag block in
`mod_portho.c` is commented out — the live distributed path is
`dist_linalg.c` / `run_dist_postfilter`.

## Key data structures

State storage is a flat `double*` of length
`n_states × nspinngrid × complex_idx` (`complex_idx = 2` when complex,
`nspinngrid = nspin·nx·ny·nz`); each state occupies a contiguous block. The main
structs (all in `fd.h`):

- **`par_st par`** — energy windows, `dt`, `dE`, `Emin`, `Emax`.
- **`index_st ist`** — all counts and sizes (`ngrid`, `nspinngrid`,
  `complex_idx`, `mn_states_tot`, `n_filters_per_rank`, …).
- **`parallel_st`** — MPI rank/size and k-group info.
- **`grid_st`** — the real and k-space grids and the volume element `dv`.
- **`zomplex`** — `{ double re, im }`.
