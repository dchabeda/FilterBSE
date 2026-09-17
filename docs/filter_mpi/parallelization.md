# Parallelization scheme and tips

`filter_mpi` is an MPI code. The filter-diagonalization pipeline has three stages —
**(1) Filter**, **(2) Orthogonalize**, **(3) Diagonalize** — and all three are
MPI-distributed. This page explains the data layout, the communication pattern,
and the practical knobs you need to pick sensible rank counts.

!!! info "Where this is implemented"
    The distributed orthogonalization and diagonalization live in
    `dist_linalg.c` / `dist_linalg.h`, wired into `main.c` behind the
    `flag.MPIOrtho` path (non-periodic runs; the periodic per-k path is
    unchanged). Despite living on a branch once named `scalapack`, the
    implementation deliberately uses **no ScaLAPACK/BLACS** — only LAPACK
    `zheev` / `cblas_zgemm` on small `N×N` matrices plus MPI collectives. This
    sidesteps ScaLAPACK linking pain on Perlmutter; the only extra object to
    link is `dist_linalg.o`.

## Data layout: column (state-block) distribution

Let the matrix of filtered states be

\[
  A = \big[\,a_0\ \ a_1\ \ \cdots\ \ a_{N-1}\,\big]\in\mathbb{C}^{M\times N},
  \qquad M=\Ngrid \sim 10^{8},\quad N \sim 10^{3},
\]

where each **column** \(a_i\in\mathbb{C}^{M}\) is one full-grid wavefunction and
\(M=\texttt{nspinngrid}\) is the (huge) number of grid degrees of freedom. The
matrix is extremely **tall and skinny**. With \(P\) MPI ranks, rank \(p\) owns a
contiguous block of \(n_p\) columns,

\[
  A_p = A[:,\,\mathcal{C}_p]\in\mathbb{C}^{M\times n_p},\qquad
  \mathcal{C}_p=[\,o_p,\,o_p+n_p\,),\qquad \sum_{p=0}^{P-1} n_p = N,
\]

so global state \(g\) lives on rank \(\lfloor g/n\rfloor\). This is exactly the
distribution the **Filter** stage produces, and exactly what the Hamiltonian
application needs — each FFT-based \(\hat H|\psi\rangle\) requires one *whole*
state on one rank. It is preserved throughout the pipeline.

!!! note "The key invariant"
    The full \(A\) (\(\sim\!1.6\) TB complex at production sizes) is **never**
    held on a single node. Rank \(p\) stores only its \(n_p\cdot M\) complex
    numbers — \(1/P\) of the full matrix. Whole states (columns) stay put; only
    tiny \(N\times N\) matrices are ever gathered.

The physical inner product carries the grid volume element \(\dv=\mathrm{d}x\,\mathrm{d}y\,\mathrm{d}z\):

\[
  \langle x,\,y\rangle = \dv\sum_{g=0}^{M-1} \overline{x_g}\,y_g,
  \qquad\text{so}\qquad
  S = \dv\,A\Herm A \in\mathbb{C}^{N\times N},\quad S_{ij}=\langle a_i,\,a_j\rangle .
\]

## The two ring primitives

Both distributed stages reduce to two operations on column-distributed data.
Each circulates the large column blocks once around a rank ring
(`MPI_Sendrecv_replace`), doing a local `zgemm` at every hop, while only tiny
\(N\times N\) data is ever reduced.

=== "(P1) Distributed Gram / overlap matrix"

    \(S=\dv\,A\Herm B\) with \(A,B\) identically column-distributed. Rank \(p\)
    holds its panel \(B_p\) **fixed** and receives every \(A_q\) in turn; the
    local product

    \[
      S[\mathcal{C}_q,\mathcal{C}_p] = \dv\,A_q\Herm B_p \in\mathbb{C}^{n_q\times n_p}
    \]

    (one `zgemm`, contraction over the \(M\) rows) fills the column panel
    \(S[:,\mathcal{C}_p]\). After \(P\) hops an `MPI_Allgatherv` of the panels
    assembles the full Hermitian \(S\) on **every** rank. For orthogonalization
    \(B=A\); for diagonalization \(B=\hat H U\).

=== "(P2) Distributed back-transform"

    \(\mathrm{out}=A\,W\) with \(W\in\mathbb{C}^{N\times r}\) small and
    replicated. Rank \(p\) accumulates its output columns by receiving each
    \(A_q\) and adding

    \[
      \mathrm{out}[:,\mathcal{K}_p] \mathrel{+}= A_q\, W[\mathcal{C}_q,\mathcal{K}_p]
    \]

    (one accumulating `zgemm` per hop). The output is left column-distributed,
    ready for the next stage.

**Communication cost.** The large column blocks circulate once (\(P-1\) neighbour
exchanges); total data moved is \(O(MN)\) — the full matrix once, *independent of
\(P\)*. Local flops per primitive are \(O(MN^2/P)\). Because \(N\sim 10^3\), the
\(O(MN^2)\) compute dominates the \(O(MN)\) communication by a factor \(\sim N\):
**the ring exchanges are cheap relative to the GEMMs.**

## Stage 2 — Distributed orthogonalization (SVD via the Gram matrix)

The goal is an orthonormal basis of the filtered states, discarding numerically
redundant directions with the singular-value threshold `SVDEPS`. For an \(M\gg N\)
matrix, forming the Gram matrix and diagonalizing it is the standard tall-skinny
SVD and maps perfectly onto the column layout.

1. **Gram matrix (P1).** Form \(S=\dv\,A\Herm A\) via the ring. Since
   \(A\Herm A = V\Lambda V\Herm\) shares its right singular vectors and squared
   singular values with the SVD \(A=U_A\Sigma V\Herm\), one small Hermitian
   eigensolve gives full spectral information: \(\sigma_i=\sqrt{\lambda_i}\).
   Solved with LAPACK `zheev` on the root and broadcast (every rank
   bit-identical).
2. **Rank / cutoff.** Retain the leading \(r\) directions with
   \(\sigma_i/\sigma_1 \ge \varepsilon\), where
   \(\varepsilon=\max(\texttt{SVDEPS},\ \varepsilon_{\mathrm{floor}})\),
   \(\varepsilon_{\mathrm{floor}}=10^{-7}\).
3. **Orthonormal basis (P2).** Build \(W=V_r\Lambda_r^{-1/2}\) and apply
   \(U=A\,W\) via the ring. The columns of \(U\) are orthonormal in the weighted
   inner product.
4. **Second pass (polish).** Repeat once with \(A\leftarrow U\) (a
   CholeskyQR2-style refinement). The second Gram is \(\approx I_r\), so it does
   not truncate further but restores \(\dv\,U\Herm U = I_r\) to machine precision
   (\(\lesssim 10^{-15}\)).

!!! warning "`DIST_SVD_FLOOR = 1e-7` — the distributed path keeps slightly fewer states"
    A Gram (normal-equations) method **squares** the condition number, so it
    cannot resolve singular values below
    \(\sqrt{\varepsilon_{\mathrm{mach}}}\approx 1.5\times10^{-8}\) relative to
    \(\sigma_1\). The effective cutoff is therefore clamped to
    \(\max(\texttt{SVDEPS}, 10^{-7})\), whereas the serial code's `SVDEPS`
    default is \(10^{-10}\). The distributed path thus keeps a few fewer states
    than serial — but the trimmed tail is exactly the ill-conditioned
    directions the \(\sigma_E\) ghost filter would remove anyway. All converged
    physical eigenvalues agree with the serial reference to \(\sim\!10^{-6}\) eV.
    Exact \(10^{-10}\) resolution at scale would require a QR-based TSQR instead
    of the Gram matrix (not built).

## Stage 3 — Distributed diagonalization

The eigenvectors are computed in the small orthonormal basis \(U\) and
transformed back to the grid. Every large operation is again local + one small
reduction.

1. **Apply the Hamiltonian (local).** Each rank applies \(\hat H\) (kinetic via
   FFT, local potential, non-local / spin–orbit projectors) to *its own* states,
   reusing `p_hamiltonian`. Embarrassingly parallel.
2. **Reduced Hamiltonian (P1).** Form the \(r\times r\) matrix
   \(\widetilde H = \dv\,U\Herm \hat H\,U\) with the Gram primitive
   (\(B=\hat H U\)).
3. **Small eigenproblem.** On the root then broadcast:
   \(\widetilde H = C\,E\,C\Herm\), giving quasiparticle energies \(E\) and
   eigenvectors \(C\) in the filtered basis.
4. **Back to the grid basis (P2).** \(\Psi = U\,C\) via a ring back-transform,
   leaving \(\Psi\) column-distributed.
5. **Eigenvalue variance (local).** The ghost-state diagnostic is purely local
   because each rank owns whole states. Using Hermiticity
   \(\langle\psi|\hat H^2|\psi\rangle = \|\hat H\psi\|^2\), so a **single** extra
   Hamiltonian apply suffices:
   \(\sigma_{E,k}^2 = \|\hat H\psi_k\|^2 - E_k^2\).

    !!! danger "One H apply, not two"
        \(\sigma_E\) uses \(\|\hat H\psi\|^2\), a *single* Hamiltonian
        application. An early bug applied \(\hat H\) twice and computed
        \(\langle\psi|\hat H^3|\psi\rangle\) instead — if you touch this code,
        keep it to one apply.

**Output.** Eigenvectors are written to `psi.dat` with **collective MPI-IO**:
state \(g\) is placed at byte offset \(g\cdot\texttt{stlen}\cdot 8\) with
\(\texttt{stlen}=\texttt{complex\_idx}\cdot M\), matching the serial writer and
the `get_n_states.x` / BSE readers. Each rank writes its block directly at its
global offset; energies and \(\sigma_E\) are gathered to the root for `eval.dat`.

## Cost summary

| Operation | Local flops | Communication | Reduced size |
|---|---|---|---|
| Gram \(S=\dv A\Herm A\) (P1) | \(O(MN^2/P)\) | \(O(MN)\) ring + Allgatherv | \(N\times N\) |
| Back-transform \(A\,W\) (P2) | \(O(MNr/P)\) | \(O(MN)\) ring | — |
| Small eigensolve (`zheev`) | \(O(N^3)\) (root) | Bcast \(N\times N\) | \(N\times N\) |
| \(\hat H U\) (per state) | \(O(\tfrac{N}{P} M\log M)\) | none | — |

Memory per rank stays at a small multiple of \(A_p\); no operation ever
materialises more than \(O(MN/P)\) (one rank's share) of the large data.

## Practical tips for choosing ranks

- **Rank count vs. states.** States are split into contiguous blocks of columns,
  so with \(P\) ranks each owns \(\approx N/P\) states. Pick \(P\) so that
  \(N/P\) states fit comfortably in a rank's memory alongside working buffers
  (a small multiple of one state block). Very large \(P\) with few states leaves
  ranks idle — there is no benefit to \(P > N\).
- **Filter is embarrassingly parallel** across the random starting vectors, so it
  scales well; ortho/diag are dominated by local GEMMs (\(O(MN^2/P)\)) with cheap
  \(O(MN)\) ring communication, so they also scale with \(P\) until \(N^3\)
  root-side `zheev` or the ring latency starts to matter (only at very large
  \(P\), small \(N\)).
- **Use the distributed path (`MPIOrtho`) when the full state matrix would not
  fit on one node.** For small toy systems the serial gather path is fine and
  slightly more accurate at the \(10^{-10}\) tail.
- **Restart at scale.** The distributed writer produces `psi.dat` + `eval.dat`
  (and `psi.par` / `eval.par` / `conf.par`); run BSE from these via its "unsafe"
  input path rather than expecting a monolithic `output.dat`.

## Validation

| Test | Metric | Result |
|---|---|---|
| Unit: ring Gram vs. dense (4 ranks) | \(\max|S_{\text{ring}}-S_{\text{ref}}|\) | \(1.6\times10^{-13}\) |
| Unit: distributed SVD, rank-deficient | retained rank \(r\) | correct (uneven split) |
| Unit: orthonormality | \(\max|\dv\,U\Herm U-I|\) | \(2.4\times10^{-15}\) |
| End-to-end CsPbI\(_3\) (4 ranks) vs. serial | converged eigenvalues | match to \(1.6\times10^{-6}\) eV |
| End-to-end: `psi.dat` | byte layout / offsets | exact |

The self-test can be compiled with `-DDIST_LINALG_TEST` in `dist_linalg.c` and
run on \(N\) ranks.
