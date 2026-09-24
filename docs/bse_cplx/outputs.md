# Outputs

`bse_cplx` writes ASCII (columnar) text for all physics results; wavefunction and
kernel I/O is raw binary. Energies are in atomic units unless a column is labeled
eV (\(\texttt{AUTOEV} = 27.2114\)). Rank 0 does all file I/O.

## Exciton spectrum

### `exciton.dat`

The exciton spectrum, one line per exciton:

```
#n   E_n   <H>   <H_dir>   <H_exc>   <H_0>   E_B(eV)
```

\(E_B\) is the binding energy \((E_n - \langle H_0\rangle)\cdot\texttt{AUTOEV}\).
The \(\langle H_{\text{dir}}\rangle\), \(\langle H_{\text{exc}}\rangle\),
\(\langle H_0\rangle\) columns are the direct, exchange, and diagonal expectation
values per exciton.

### `BSEcoeff.dat`

The full set of exciton coefficients (BSE eigenvectors), written by both the
serial and distributed (ScaLAPACK) solve paths. One row per two-particle basis
component \(i\), one column per exciton \(j\); each entry is the complex
coefficient \(C_{i}^{j}\) written as a `re, im` pair (row `i`, column `j` holds
component `i` of exciton `j`). Consumed by `filter/carrier_dens.c`
(`fscanf("%lg, %lg")`).

### `BSEeval.par`

The selected QP basis (from stage 1): `idx  eigval  sigma` per state.

### Fixed-point exciton densities (`rhoElec<n>.cube`, `rhoHole<n>.cube`)

Set `nExcitonDensities = M` to write Gaussian cube files of the fixed-point
carrier densities for the lowest `M` excitons (`n = 0 … M-1`). For each exciton
two cubes are written:

- `rhoElec<n>.cube` — **electron** carrier density (hole integrated out):
  \(\rho^e_n(\mathbf r_e) = \sum_i \sum_\sigma \big|\sum_a C^n_{ia}\,\phi^\sigma_a(\mathbf r_e)\big|^2\)
- `rhoHole<n>.cube` — **hole** carrier density (electron integrated out):
  \(\rho^h_n(\mathbf r_h) = \sum_a \sum_\sigma \big|\sum_i C^n_{ia}\,\phi^\sigma_i(\mathbf r_h)\big|^2\)

Integrating out the other carrier collapses the two-particle sum by QP
orthonormality (\(\int d\mathbf r\,\phi_i^*\phi_j = \delta_{ij}\)), so each density
is real and non-negative. The cube value is the **signed density** rendered with
the same method as `filter_mpi/makecube`
(\(\sqrt{\mathrm{Re}^2+\mathrm{Im}^2}\cdot\mathrm{sign}\)), not the squared
density. Written by rank 0 as a post-processing step after the BSE solve (needs
`conf.par` in the run directory, the same config file the rest of the run
reads); works in both the serial and distributed paths.

## Optical properties (exciton)

| File | Contents |
|---|---|
| `OS.dat` | electric oscillator strengths: `n  sqrt(|μ|²)  E_n  f_osc  Re/Im(μx,μy,μz)` |
| `M.dat` | magnetic dipole / magnetic OS (same layout, \((4/3)E|m|^2\)) |
| `rs.dat` | exciton rotational strengths (CD): `n  E_n  R_n`, with the full bilinear \(R_n=\mathrm{Im}[\boldsymbol\mu_n\cdot\mathbf m_n^*]\) including electron–hole interference; conserves \(\sum_n R_n=\sum_{ia} r_{ia}\) (sum-rule line printed to stdout) |
| `bs-coeff-<n>.dat` | per-exciton BSE coefficients (for the first 10 excitons): `hole  elec  E_ia  Re(C)  Im(C)  |C|²` |

## Single-particle (QP) properties

| File | Contents |
|---|---|
| `OS0.dat` | single-particle electric transition dipoles: `i  a  sqrt(mu2)  Ea-Ei  f_osc  μ components…` |
| `M0.dat` | single-particle magnetic dipoles |
| `rs0.dat` | single-particle rotational strengths |
| `qp_spins.dat` | per-QP-state spin-up/down fractions (spinor runs) |

## Angular momentum (`calcSpinAngStat = 1`)

**QP-basis matrices** (each line `i  j  Re  Im`, holes then electrons):

- `sx.dat`, `sy.dat`, `sz.dat` — spin matrices
- `ssqr.dat` — \(S^2\) (diagonal = \(3/4\) for spin-½ QP states; a normalization/orthonormality check)
- `lx.dat`, `ly.dat`, `lz.dat` — orbital angular-momentum matrices
- `lsqr.dat` — \(L^2\)
- `ls.dat` — \(L\cdot S\)

**Exciton-basis:**

- `spins.dat` — `n  Re/Im(Sx,Sy,Sz)  S_tot`
- `orbital.dat` — `n  Lx  Ly  Lz  L_tot`
- `couple.dat` — `n  <L·S>`

The QP-basis grid integrals (both the spin matrices and the angular-momentum
matrix elements) are offloaded to the GPU when `gpuAccel = 1` and a device is
visible, using the same resident-state technique as the Coulomb kernel (the
FFT-based \(L|\psi\rangle\) operators stay on the CPU); the result is
bit-identical to the CPU path and it falls back to CPU automatically.

## Kernel checkpoints

| File | Contents |
|---|---|
| `direct-<rank>.dat`, `exchange-<rank>.dat` | per-rank ASCII kernel dumps (`a b i j ibs jbs Re Im`), used by `restartCoulomb` |
| `direct.dat`, `exchange.dat` | whole-matrix kernel reload, used by `coulombDone = 1` |

See [Job options → Restarting](job-options.md#restarting-and-checkpointing).

## Debug matrix dumps

- `h0.dat` — the diagonal \(h_0\) matrix
- `bsRE.dat`, `bsIM.dat` — the assembled BSE matrix (real/imag)
- `HBSmatRE.dat`, `HBSmatIM.dat` — BSE matrix (serial path only)

## stdout / stderr

Rank 0 prints a banner and timestamped stage headers ("1. INITIALIZING JOB",
"2. COMPUTING … POTENTIALS", "3. COMPUTING … KERNEL", "4. …"), plus:

- QP basis diagnostics: HOMO/LUMO indices, energies, fundamental gap
- the GPU vs CPU path chosen, and the block-cyclic layout summary
- the ground-state exciton energy
- kernel and diagonalization wall times, and a final CPU/wall-time summary

!!! note
    The stage numbers on stdout are slightly inconsistent — both the dipole and
    BSE stages print "4." This is cosmetic. Errors go to stderr; the example
    submit script routes `2> error.dat > run.dat`.
