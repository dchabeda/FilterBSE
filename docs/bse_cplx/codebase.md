# Codebase & module layout

`bse_cplx` is C + MPI + OpenMP with **native `double complex`** typing throughout.
It links either cray-libsci (ScaLAPACK + LAPACKE, the default GPU/distributed
build) or Intel MKL (the CPU-only build), plus FFTW. See
[Job options → Building](job-options.md#building) for the two Makefiles and their
executables (`bse_cplx_gpu.x`, `bse_cplx.x`).

## Top-level control flow (`main.c`)

`main()` runs the pipeline as a sequence of module drivers:

1. `MPI_Init`; build a **node-local shared communicator** via
   `MPI_Comm_split_type(MPI_COMM_TYPE_SHARED)`.
2. Read `argv[1]` → `initUnsafe` (safe vs unsafe input path).
3. **`mod_init`** — read input, select the QP basis window, allocate and populate
   the node-shared `psi_qp` wavefunctions.
4. If spinors: `qp_spin_frac` (per-state spin fractions) on rank 0.
5. **`mod_dipole`** — single-particle electric/magnetic dipoles (collective).
6. If `noCalcExciton`: print and exit.
7. **`mod_pot`** — build the Coulomb potentials on the grid.
8. **`mod_kernel`** — build the direct/exchange e–h kernel.
9. **`mod_bse`** — assemble and diagonalize the BSE Hamiltonian; free the kernel.
10. On rank 0: **`calc_optical_exc`** — exciton oscillator/magnetic/rotational
    strengths.
11. If `calcSpinAngStat`: QP and exciton spin / angular-momentum matrices
    (collective; `bs_coeff` broadcast first on the distributed path).
12. Teardown: release the shared window (`MPI_Win_unlock_all` / `MPI_Win_free`),
    free arrays, `MPI_Finalize`.

!!! note "Cosmetic: two 'stage 4' labels"
    Both `mod_dipole` and `mod_bse` print a "4." stage header on stdout — a
    known cosmetic inconsistency, not a functional issue.

## Module drivers (`mod_*`)

| File | Role |
|---|---|
| `mod_init.c` | read filter output + `input.par`, select QP window, allocate node-shared `psi_qp`, load wavefunctions |
| `mod_pot.c` | allocate `pot_bare` / `pot_screened`, call `init_elec_hole_kernel` |
| `mod_kernel.c` | allocate direct/exchange (full N² or block-cyclic tile), dispatch `calc_eh_kernel_cplx` or load from disk; handle `calcCoulombOnly` |
| `mod_bse.c` | assemble \(H = h_0 - \text{direct} - \text{exchange}\); `build_h0_mat`, `build_BSE_mat` (serial path) |
| `mod_dipole.c` | drive the single-particle dipole calculations |

## Compute & support files

| File | Role |
|---|---|
| `read.c` | parse `input.par` (`read_input`) and `unsafe_input.par` (`read_unsafe_input`); typed key parser `read_field` |
| `save.c` | `read_filter_output` (parse monolithic `output.dat`), input-state echo, density/cube printers |
| `basis.c` | `get_qp_basis_indices` (QP window selection) and `get_qp_basis` (load selected wavefunctions) |
| `init.c` | `init_elec_hole_kernel` — build bare & screened Coulomb on the reciprocal grid; `calc_coulomb` = erf(γr)/r |
| `coulomb.c` | `calc_eh_kernel_cplx` — the e–h kernel: direct + exchange, MPI parity split, **GPU offload + auto-tiling**, block-cyclic routing, restart |
| `hartree.c` | `hartree` — FFT convolution of a density with a reciprocal-space potential |
| `pbse.c` | distributed ScaLAPACK path: BLACS/descriptor setup, block-cyclic router, `bethe_salpeter_dist` (`pzheevd`/`pzhemm`/`pzgemr2d`) |
| `bethe-salpeter.c` | `bethe_salpeter` — serial single-rank assemble → `diag` → expectation values → `exciton.dat` |
| `diag.c` | `diag` — dense host eigensolve via `LAPACKE_zheev` (row-major, all vectors) |
| `dipole.c` | `calc_elec_dipole`, `calc_mag_dipole`, `calc_rotational_strength` (single-particle) |
| `angular.c` | QP and exciton spin / orbital angular-momentum matrices (`calc_qp_*`, `calc_xton_*`) |
| `ang_mom_operators.c` | `l_operator` (angular momentum via FFT k-space derivatives), `p_operator`, G-vector/k-vector setup |
| `spin.c` | `qp_spin_frac` — per-state spin-up/down fraction |
| `optical.c` | `calc_optical_exc` — exciton oscillator strengths, magnetic OS, rotational strength, per-exciton coefficients |
| `norm.c` | normalization helpers |
| `aux.c` | timing, aligned allocation (the `ALLOCATE` macro), dynamic workload chunking, progress bars |
| `write.c` | cube files, stdout banners, state dumps |
| `nerror.c` | `nerror`, `terminate` fatal-error helpers |

## Headers

- **`fd.h`** — the master header: all structs (`index_st ist`, `grid_st`,
  `par_st`, `flag_st`, `parallel_st`), macros (`ALLOCATE`, `AUTOEV = 27.2114`,
  `PR_LEN`), and the prototypes for modules without their own header
  (`read_filter_output`, `calc_optical_exc`, `diag`, …). It includes `<lapacke.h>`
  under `USE_LIBSCI` or `<mkl.h>` otherwise.
- Per-module headers: `main.h`, `mod_init.h`, `mod_pot.h`, `mod_kernel.h`,
  `mod_bse.h`, `mod_dipole.h`, `read.h`, `basis.h`, `init.h`, `coulomb.h`,
  `hartree.h`, `bethe-salpeter.h`, `dipole.h`, `angular.h`, `spin.h`, `aux.h`,
  `write.h`, `pbse.h`. (There is no `save.h` / `optical.h` / `diag.h` — those
  prototypes live in `fd.h`.)

## Not in the active build

- `coulomb_same.c` — a duplicate/older copy of `calc_eh_kernel_cplx`; **not** in
  any Makefile's object list. Do not treat it as live.
- `convert.c`, `merge_psi_spinor.c` — standalone tools with their own `main`, not
  part of the solver.
- `test.c` — empty.
- `pbse.c` also contains a `-DPBSE_TEST` standalone `main` that self-tests the
  distributed eigensolve against serial `LAPACKE_zheev`.

## Key data structures

| Name | What it holds |
|---|---|
| `psi_qp` | selected QP wavefunctions, `double complex[n_qp × nspinngrid]`; in an MPI-3 shared window (one copy per node) |
| `eig_vals`, `sigma_E` | QP energies and variances (from the filter output) |
| `pot_bare`, `pot_screened` | reciprocal-space Coulomb potentials, `double complex[ngrid]` |
| `direct`, `exchange` | the e–h kernel matrices; full `n_xton²` (serial) or a block-cyclic tile (distributed) |
| `xton_ene` | exciton energies (all ranks) |
| `bs_coeff` | exciton eigenvectors, `double complex[n_xton²]`; component *i* of exciton *j* at `bs_coeff[i·n_xton + j]` |

`n_xton = n_holes × n_elecs` is the BSE matrix dimension.
