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

### `BSEeval.par`

The selected QP basis (from stage 1): `idx  eigval  sigma` per state.

## Optical properties (exciton)

| File | Contents |
|---|---|
| `OS.dat` | electric oscillator strengths: `n  sqrt(|μ|²)  E_n  f_osc  Re/Im(μx,μy,μz)` |
| `M.dat` | magnetic dipole / magnetic OS (same layout, \((4/3)E|m|^2\)) |
| `rs.dat` | rotational strengths: `n  E_n  R` |
| `bs-coeff-<n>.dat` | per-exciton BSE coefficients (for the first 10 excitons): `hole  elec  E_ia  Re(C)  Im(C)  |C|²` |

## Single-particle (QP) properties

| File | Contents |
|---|---|
| `OS0.dat` | single-particle electric transition dipoles: `i  a  sqrt(mu2)  Ea-Ei  f_osc  μ components…` |
| `M0.dat` | single-particle magnetic dipoles |
| `rs0.dat` | single-particle rotational strengths |
| `qp_spins.dat` | per-QP-state spin-up/down fractions (spinor runs) |

## Angular momentum (`calcSpinAngStat = 1`)

**QP-basis matrices:**

- `sx.dat`, `sy.dat`, `sz.dat` — spin matrices
- `lx.dat`, `ly.dat`, `lz.dat` — orbital angular-momentum matrices
- `lsqr.dat` — \(L^2\)
- `ls.dat` — \(L\cdot S\)

**Exciton-basis:**

- `spins.dat` — `n  Re/Im(Sx,Sy,Sz)  S_tot`
- `orbital.dat` — `n  Lx  Ly  Lz  L_tot`
- `couple.dat` — `n  <L·S>`

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
- `BSEcoeffRE.dat`, `BSEcoeffIM.dat` — exciton coefficients (serial path only)

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
