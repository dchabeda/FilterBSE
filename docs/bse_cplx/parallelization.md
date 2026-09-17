# Parallelization scheme and tips

`bse_cplx` combines four kinds of parallelism: an **MPI parity split** of the
Coulomb kernel, **ScaLAPACK** for the distributed dense BSE eigensolve, an
**MPI-3 node-shared** copy of the wavefunctions, and optional **GPU offload** of
the kernel with automatic memory tiling. All of it is active only in the
`bse_cplx_gpu.x` build (`-DUSE_SCALAPACK -DUSE_GPU_OFFLOAD`); the MKL build
(`bse_cplx.x`) is single-node serial.

!!! info "Automatic dispatch"
    On the distributed build, `dist = (mpi_size > 1)`. With more than one rank you
    get the block-cyclic kernel + `pzheevd`; with one rank (or the MKL build) you
    get the full-matrix kernel + serial `LAPACKE_zheev`. You do not toggle this
    yourself.

## MPI decomposition of the Coulomb kernel

The kernel (stage 2) is the expensive part, so it is distributed two ways at once:

- **Parity split (direct vs exchange).** Ranks are split by
  `rank % 2` into two sub-communicators: **even ranks compute the direct
  (screened) term \(K^d\); odd ranks compute the exchange (bare) term \(K^x\)** —
  concurrently. This is why you want an **even number of ranks ≥ 2**; a single
  rank does not exercise the split correctly.
- **Work distribution within a group.** Inside each parity group the
  electron-pair index space is distributed strided over the ranks
  (`start = rank_in_group; step = group_size`), and the inner hole loops are
  OpenMP-threaded.

The code exploits the two-fold permutation symmetry of the complex two-electron
integrals (\([ij|ab] = [ji|ba]^*\)), so only the lower triangle is computed and
stored.

## Kernel storage: block-cyclic vs full matrix

At large basis sizes the kernel matrices dominate memory, so how they are stored
matters:

=== "Distributed (mpi_size > 1)"

    Each rank allocates only its **block-cyclic ScaLAPACK tile**
    (\(\sim n_{\text{xton}}^2 / P\)) of `direct` and `exchange`. Computed elements
    are binned by block-cyclic owner into a **router** and delivered with a single
    `MPI_Alltoallv`. **No rank ever holds the full matrix** — adding ranks
    genuinely reduces per-rank memory, and the tile is consumed in place by
    `pzheevd`.

    !!! note "Why this exists"
        Previously `direct` and `exchange` were each allocated as the full
        \(n_{\text{xton}}^2\) on *every* rank (e.g. 16.8 GB each at
        \(n_{\text{xton}} = 32400\)) and then reduced to rank 0 — so adding nodes
        never helped. The block-cyclic layout fixed that OOM.

=== "Serial (single rank / MKL build)"

    The full \(n_{\text{xton}}^2\) matrices are allocated; even/odd ranks
    `MPI_Reduce` their contributions to rank 0/1 and exchange is sent to rank 0.

## The BSE eigensolve: ScaLAPACK

The dense Hermitian BSE solve (stage 3) is distributed with **ScaLAPACK** on the
default build:

- A near-square **BLACS** 2D process grid, block size 64, with block-cyclic
  descriptors (`descinit_`, `numroc_`).
- **`pzheevd`** (divide-and-conquer Hermitian eigensolver, `jobz='V'`,
  `uplo='L'`) computes all eigenvalues and eigenvectors. Eigenvalues are
  broadcast to all ranks; eigenvectors are gathered to rank 0 (`pzgemr2d`) as
  `bs_coeff` for the optical/angular analysis.
- Expectation values \(\langle H_{\text{dir}}\rangle\),
  \(\langle H_{\text{exc}}\rangle\), \(\langle H_0\rangle\) are computed with
  `pzhemm` + a local contraction.
- The libsci ScaLAPACK/PBLAS are internally OpenMP-threaded (the `_mp` variant),
  so each rank uses its full `nThreads` quota.

!!! warning "This is a CPU eigensolve"
    ScaLAPACK `pzheevd` runs across **all MPI ranks' CPUs** — which is exactly the
    point (use all the cores), but it is not a GPU eigensolver. A true GPU dense
    eigensolver (cuSOLVERMp / ELPA-GPU) is not part of libsci and is future work.
    GPU offload here targets the **Coulomb kernel only**; the eigensolve stays on
    the host in every build (see [why the full spectrum is
    needed](pipeline.md#stage-3-bse-eigensolve)).

    32-bit ScaLAPACK integers are used, which is safe because only 2D indices (not
    the flat \(N^2\) count) are ever passed — good for \(n_{\text{xton}}\) up to
    ~40000.

## Distributed dipole and angular momentum

Stage 4 is also distributed (following the same idiom): each rank computes a
disjoint round-robin slice of the outer loop
(`for x = rank; x < n; x += mpi_size`), threaded with OpenMP over the grid, then
`MPI_Allreduce(SUM)` reassembles the result; rank 0 does the file I/O. This
applies to `calc_elec_dipole` (over holes), `calc_mag_dipole` (over electrons),
and the QP and exciton spin/angular routines in `angular.c`. Because every matrix
element has a single owning rank, the Allreduce reassembles disjoint pieces with
no reordering — the QP-level results are **bit-identical for any rank count**.

The exciton-basis routines read the full `bs_coeff`, which is broadcast row-by-row
to all ranks first (via an MPI datatype so the element count stays within 32-bit
range at large \(N\)).

## Node-shared `psi_qp`

The QP wavefunctions `psi_qp` are large (~0.5 GB per state). Storing one copy per
MPI rank OOMs the host at large basis sizes (e.g. 128 states × 4 ranks/node on a
256 GB node). Instead, `psi_qp` lives in an **MPI-3 shared-memory window**
(`MPI_Win_allocate_shared` over the node communicator): **one physical copy per
node**, read-only shared by all ranks on that node. Only node-rank 0 fills it from
disk; every rank resolves the same base pointer. Each rank still "sees" a full
`psi_qp` and maps its own tiles to its own GPU.

!!! note "Limit"
    Node-sharing helps until a single copy exceeds one node's RAM (roughly 256+
    states, or a very large grid). Beyond that, streaming state tiles from disk or
    distributing ownership is required — that is future work.

## GPU offload with auto-tiling

GPU offload lives entirely in `coulomb.c` (guarded by `-DUSE_GPU_OFFLOAD`), using
OpenMP `target` regions on the native `double complex` data. It is enabled by
`gpuAccel = 1` and gated at runtime by `omp_get_num_devices() > 0`.

- State residency is auto-sized from **`BSE_GPU_MEM_GB`** (default 36). If the
  required states fit, the kernel takes a fast single-map path; otherwise it
  **automatically tiles** the states (or falls back to CPU for that kernel) — no
  correctness impact.
- The direct kernel needs all holes resident; the exchange kernel needs holes and
  electrons, so it tiles more readily. Tile shapes are chosen to minimize the
  number of passes (each pass recomputes the host Hartree FFT).

!!! danger "NVHPC device-pointer pitfall (for maintainers)"
    Mapping two **disjoint slices of the same base array** into one target region
    makes NVHPC translate only the first slice → `CUDA_ERROR_ILLEGAL_ADDRESS`.
    The code copies each tile into its **own distinct contiguous buffer**
    (`drow` / `dcol`) and maps those. Distinct base pointers are fine; slices of
    one array are not. Keep it that way.

## Practical tips for choosing ranks and nodes

- **Use one MPI rank per GPU**, an even number ≥ 2, e.g.
  `srun -n 4 --gpus-per-task=1 --gpu-bind=closest -c 32`. Even count keeps the
  direct/exchange parity split balanced.
- **Add ranks to fit large kernels.** On the distributed build, per-rank kernel
  memory is \(\sim n_{\text{xton}}^2 / P\), so more ranks/nodes directly relieve
  kernel OOM.
- **Set `BSE_GPU_MEM_GB`** to `GPU_GB × fraction / ranks_per_gpu` (36 for a 40 GB
  A100, one rank/GPU). Larger budgets avoid exchange-kernel tiling and its
  redundant Hartree recomputes.
- **Set `nThreads`** to the cores available per rank; it feeds FFTW, the
  dipole/angular loops, and the threaded ScaLAPACK solve.
- **CPU-only run:** use the MKL build (`bse_cplx.x`) on a single node for small
  systems, or the distributed build with `gpuAccel = 0` / `OMP_TARGET_OFFLOAD=DISABLED`.

## Validation notes

- QP-level outputs (`OS0.dat`, `M0.dat`, the `sx/lx/…` matrices) depend only on
  `psi_qp` and are **bit-identical across any rank count** — the cleanest
  validation anchors.
- Exciton-level outputs depend on `bs_coeff`, whose per-vector phase and
  degenerate-subspace rotation are gauge-dependent under `pzheevd`. Compare
  gauge-invariant quantities (eigenvalues, \(|C|^2\), oscillator strengths,
  diagonal \(\langle n|O|n\rangle\)) to ~\(10^{-10}\), not bit-for-bit.
- The GPU and CPU Coulomb kernels are bit-identical: run once with `gpuAccel=1`
  and once with `OMP_TARGET_OFFLOAD=DISABLED` and diff `BSEeval.par`.
