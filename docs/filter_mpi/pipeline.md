# Pipeline: filter → ortho → diag

`filter_mpi` computes quasiparticle states in three stages. This page explains
what each stage does and how data flows between them. For the **distributed
(MPIOrtho)** implementation of stages 2 and 3 — the ring/Gram scheme — see
[Parallelization](parallelization.md); this page describes the pipeline
conceptually and the serial path.

```
   input.par / conf.par / pots
              │
              ▼
   ┌───────────────────────┐
   │  1. FILTER            │  Chebyshev/Newton filter of random vectors
   │  mod_filter.c         │  → redundant set of filtered states
   │  filter.c, coeff.c    │
   └──────────┬────────────┘
              │  psi_rank  (many, non-orthogonal)
              ▼
   ┌───────────────────────┐
   │  2. ORTHOGONALIZE     │  SVD; drop redundant directions
   │  mod_ortho.c/ortho.c  │  → orthonormal basis (fewer states)
   │  (or dist_linalg.c)   │
   └──────────┬────────────┘
              │  U  (orthonormal)
              ▼
   ┌───────────────────────┐
   │  3. DIAGONALIZE       │  build H in the basis, diagonalize
   │  mod_diag.c/Hmat.c    │  → energies + grid-basis eigenvectors
   │  (or dist_linalg.c)   │  + sigma_E ghost-state check
   └──────────┬────────────┘
              │
              ▼
       eval.dat / psi.dat / output.dat
```

## Stage 1 — Filter

**Files:** `mod_filter.c`, `filter.c`, `coeff.c`, `energy.c`.

The filter projects random starting vectors onto a set of energy targets inside
the valence and conduction windows, amplifying the components near each target.

1. **Energy range.** `get_energy_range` (`get_energy_range_k` for periodic)
   computes the spectral bounds `par->Emin/Emax`. For periodic systems this is
   k-dependent because the kinetic term carries \(k\cdot G\).
2. **Filter coefficients.** `gen_newton_coeff` (`coeff.c`) builds Newton
   interpolation coefficients `an` at Chebyshev support points `zn` that
   approximate a near-delta filter function at each of the `mStatesPerFilter`
   energy targets. `readCoeffs = 1` reads them from disk instead.
3. **Initial states.** `init_filter_states` fills `psi_rank` with random initial
   states — one per filter cycle, seeded via `Randomize()` or a fixed
   `rand_seed`, offset by rank.
4. **Filter cycles.** `run_filter_cycles` → `filter_cycle` runs the core
   Chebyshev/Newton recursion. Term 0 seeds each target with `an[0]·ψ`; then for
   each Chebyshev term the code applies the scaled, normalized Hamiltonian once
   (`p_hnorm`, which wraps the FFT-based kinetic + local + NL + SO Hamiltonian and
   rescales the spectrum into \([-1,1]\)) and accumulates `an[...]·ψ` into **every
   target state simultaneously**. So one Hamiltonian application serves all
   `mStatesPerFilter` targets.

**Output of stage 1:** each rank holds
`n_states_per_rank = n_filters_per_rank × mStatesPerFilter` filtered states in
`psi_rank`. These are **not orthogonal** and are redundant (many filter cycles
converge to the same physical states). Energies are written to `ene-filt-*.dat`.

!!! info "This stage is embarrassingly parallel"
    Filter cycles are independent, so each rank filters its own share
    (`nFilterCycles / mpi_size`) with no communication. See
    [Parallelization](parallelization.md).

## Stage 2 — Orthogonalize

**Files:** `mod_ortho.c` / `ortho.c` (serial) or `dist_linalg.c` (distributed).

The redundant filtered states are reduced to an orthonormal basis via an SVD.

- **Serial path.** After the filtered states are gathered to rank 0 (`psitot`),
  optional time reversal is applied, then `ortho_real` (`dgesvd`) or `ortho_cplx`
  (`zgesvd`) takes the SVD of the `ngrid × mn_states` state matrix. The number of
  retained vectors, `cutoff`, is the first index where
  \(\sigma_i/\sigma_0 < \texttt{SVDEPS}\) (\(10^{-10}\)); `mn_states_tot` is
  updated to `cutoff` and `normalize_all` follows.
- **Distributed path (`MPIOrtho = 1`).** A Gram-matrix SVD done entirely across
  ranks without gathering the states, with a coarser \(10^{-7}\) cutoff. See
  [Parallelization → Stage 2](parallelization.md#stage-2-distributed-orthogonalization-svd-via-the-gram-matrix).

**Data flow:** redundant filtered states in → orthonormal basis (same grid
layout, **fewer columns**) out.

## Stage 3 — Diagonalize

**Files:** `mod_diag.c` / `Hmat.c` (serial), `Hmat_mpi.c` or `dist_linalg.c`
(distributed).

1. **Build the subspace Hamiltonian.** `diag_H` (`Hmat.c`) forms
   \(\widetilde H_{ij} = \langle\psi_i|\hat H|\psi_j\rangle\,\mathrm{d}v\) over the
   orthonormal states, applying the FFT Hamiltonian to each.
2. **Diagonalize.** `dsyev` (real) or `zheev` (complex) diagonalizes the small
   \(\widetilde H\). Eigenvalues (ascending) go into `eig_vals`; eigenvectors are
   back-transformed to the grid basis. The periodic path calls `diag_H` per
   k-point; `MPIDiag = 1` uses `diag_H_mpi`; the `MPIOrtho` path uses
   `dist_diag`.
3. **Eigenvalue variance.** `mod_sigma` / `calc_sigma_E` computes
   \(\sigma_E = \sqrt{\langle\psi|\hat H^2|\psi\rangle - \langle\psi|\hat H|\psi\rangle^2}\)
   for each state — a ghost-state / convergence diagnostic. States with
   \(\sigma_E \le\) `sigmaECut` are "converged". (In the distributed path this
   uses \(\|\hat H\psi\|^2\), a single extra Hamiltonian apply.)

**Data flow:** orthonormal basis in → grid-basis eigenvectors (`psitot` /
`psi_rank`) + `eig_vals` + `sigma_E` out → `mod_output`.

## What comes out

The final stage (`mod_output`) writes `eval.dat` (energies + \(\sigma_E\)),
`psi.dat` (grid-basis eigenvectors), and, unless on the distributed path, the
monolithic `output.dat` consumed by the BSE codes. See [Outputs](outputs.md).
