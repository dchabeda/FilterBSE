# Pipeline

`bse_cplx` turns the filter's quasiparticle (QP) states into excitons and optical
properties in four stages. This page describes each stage and how data flows
between them. For how stages 2 and 3 are distributed across MPI ranks and GPUs,
see [Parallelization](parallelization.md).

```
   filter output (output.dat  or  psi.par/eval.par/conf.par)
              │
              ▼
   ┌───────────────────────────┐
   │  1. QP BASIS              │  read + apply sigma_E / Fermi window
   │  mod_init, basis.c        │  → n_holes, n_elecs; load psi_qp
   └──────────┬────────────────┘
              │  psi_qp, eig_vals   (n_xton = n_holes × n_elecs)
              ▼
   ┌───────────────────────────┐
   │  2. e–h KERNEL           │  FFT Coulomb convolutions
   │  mod_pot, mod_kernel      │  → direct (screened) + exchange (bare)
   │  coulomb.c, hartree.c     │
   └──────────┬────────────────┘
              │  direct, exchange
              ▼
   ┌───────────────────────────┐
   │  3. BSE EIGENSOLVE       │  H = h0 − (direct + exchange)
   │  mod_bse, pbse.c/diag.c   │  → xton_ene + ALL eigenvectors bs_coeff
   └──────────┬────────────────┘
              │  xton_ene, bs_coeff
              ▼
   ┌───────────────────────────┐
   │  4. OPTICAL & ANGULAR    │  contract bs_coeff with dipoles / L, S
   │  optical.c, dipole.c,     │  → spectra + state character
   │  angular.c                │
   └──────────┬────────────────┘
              ▼
      exciton.dat / OS.dat / M.dat / spins.dat / …
```

## Stage 1 — QP basis and energy window

**Files:** `mod_init.c`, `save.c` / `read.c`, `basis.c`.

1. **Read the filter output.** `read_filter_output` (safe path) or
   `read_unsafe_input` (unsafe path) loads the grid, atoms, `eig_vals`, and
   `sigma_E`; `read_input` supplies `fermiEnergy`, `sigmaECut`, and the state
   caps. See [Input files](input-files.md).
2. **Select the window.** `get_qp_basis_indices` classifies each filter state:

    - **hole** if `sigma_E < sigmaECut && eigval < fermiEnergy`
    - **electron** if `sigma_E < sigmaECut && eigval > fermiEnergy`

    The \(\sigma_E\) cut discards unconverged/ghost states (which is why a simple
    contiguous window is wrong — see the
    [unsafe-init note](input-files.md#two-input-paths)). It then optionally keeps
    only `maxHoleStates` / `maxElecStates` nearest the gap, and records
    `homo_idx`, `lumo_idx`, `n_holes`, `n_elecs`. `BSEeval.par` is written with
    the selected basis.
3. **Load wavefunctions.** `get_qp_basis` seeks and reads only the selected
   states into `psi_qp` — stored in a node-shared window (see
   [Parallelization](parallelization.md#node-shared-psi_qp)).

**Output:** `psi_qp`, `eig_vals`, and the counts. The BSE matrix dimension is
`n_xton = n_holes × n_elecs`.

## Stage 2 — Electron–hole kernel

**Files:** `mod_pot.c` / `init.c`, `mod_kernel.c` / `coulomb.c`, `hartree.c`.

1. **Build the Coulomb potentials.** `init_elec_hole_kernel` builds two
   reciprocal-space potentials on the grid: a **bare** Coulomb (unscreened, used
   by the exchange term) and a **screened** Coulomb (dielectric-screened,
   anisotropic via `epsX/Y/Z`, Rohlfing–Louie form, used by the direct term). The
   \(1/r\) singularity is regularized as \(\operatorname{erf}(\gamma r)/r\).
2. **Form the kernel.** `calc_eh_kernel_cplx` computes, for each pair of
   electron–hole transitions, a two-particle Coulomb integral: build the density
   \(\rho = \psi_a^* \psi_b\), convolve it with a potential via an FFT
   (`hartree`), and integrate against \(\psi_i^* \psi_j\). The **direct** term
   \(K^d\) uses the screened potential; the **exchange** term \(K^x\) uses the
   bare potential.

**Output:** the `direct` and `exchange` matrices — stored either as full
`n_xton²` (serial) or as per-rank block-cyclic tiles (distributed). This stage is
the most expensive and is the target of MPI distribution and GPU offload.

## Stage 3 — BSE eigensolve

**Files:** `mod_bse.c`, `pbse.c` (distributed) / `bethe-salpeter.c` + `diag.c`
(serial).

The BSE Hamiltonian is

\[
  H = h_0 - (K^d + K^x),
\]

where \(h_0\) is diagonal with the QP transition energies \(E_a - E_i\). It is
Hermitian and dense.

- **Distributed:** `bethe_salpeter_dist` assembles \(H\) in each rank's
  block-cyclic tile and diagonalizes in place with ScaLAPACK `pzheevd`
  (`uplo='L'`), broadcasting eigenvalues and gathering eigenvectors to rank 0.
- **Serial:** `build_BSE_mat` / `build_h0_mat` then `diag` → `LAPACKE_zheev`
  (row-major, `uplo='U'`).

!!! important "The full eigenspectrum is required — this is physics, not a choice"
    `bse_cplx` computes **all** eigenvectors (`jobz = 'V'`). Optical matrix
    elements and angular-momentum expectation values are sums over the *entire*
    exciton spectrum: every exciton mixes all electron–hole pairs through its
    `bs_coeff` column, and stage 4 contracts the full coefficient set for every
    exciton. A dense eigensolver is therefore the correct algorithm — do **not**
    substitute an iterative/subset solver (Lanczos, LOBPCG, SLEPc) or drop
    eigenvectors. Distributed dense solvers (ScaLAPACK/ELPA) only become necessary
    above ~10⁵; typical runs are ~200 holes × 200 electrons ⇒ a 40000×40000
    matrix, routine on CPU.

**Output:** `xton_ene` (exciton energies, on all ranks) and `bs_coeff` (exciton
eigenvectors, gathered to rank 0). `exciton.dat` records each exciton's energy,
the \(\langle H\rangle\), \(\langle H_{\text{dir}}\rangle\),
\(\langle H_{\text{exc}}\rangle\), \(\langle H_0\rangle\) expectation values, and
the binding energy.

## Stage 4 — Optical and angular-momentum properties

**Files:** `optical.c`, `dipole.c`, `angular.c`, `spin.c`,
`ang_mom_operators.c`.

- **Single-particle dipoles** (computed earlier in `mod_dipole`):
  `calc_elec_dipole` (electric) and `calc_mag_dipole` (magnetic, via the FFT-based
  \(L\) operator) build transition dipoles between QP states.
- **Exciton optics.** `calc_optical_exc` contracts `bs_coeff` with the
  single-particle dipoles to get each exciton's dipole, oscillator strength
  \(\tfrac{2}{3} E |\mu|^2\), magnetic oscillator strength, and rotational
  strength.
- **Angular momentum** (when `calcSpinAngStat = 1`). `calc_qp_spin_mtrx` /
  `calc_qp_ang_mom_mtrx` build QP-basis \(\langle S\rangle\), \(\langle L\rangle\),
  \(\langle L^2\rangle\), \(\langle L\cdot S\rangle\); `calc_xton_spin_mtrx` /
  `calc_xton_ang_mom_mtrx` transform these into the exciton basis via `bs_coeff`.
  Spinor spin fractions come from `qp_spin_frac`.

**Output:** the spectra and state-character files described in
[Outputs](outputs.md).
