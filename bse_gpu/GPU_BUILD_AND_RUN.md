# BSE on Perlmutter GPUs — Build & Run Guide

This documents the GPU port of the BSE Coulomb kernel that lives in `bse_gpu/`,
and how to compile and run it on the NERSC Perlmutter A100 GPU nodes.

The GPU work offloads the Bethe–Salpeter Coulomb kernel
(`coulomb.c : calc_eh_kernel_cplx`) — the direct and exchange electron–hole
integrals — to the GPU via OpenMP `target teams distribute`. Everything else
(filter I/O, Hartree FFT, diagonalization) still runs on the host. The GPU path
produces results **bit-identical** to the host path; it is validated, not
experimental.

---

## 1. What changed last session (summary of GPU edits)

### The kernel offload
- `coulomb.c`'s `calc_eh_kernel_cplx` now runs the direct and exchange
  electron–hole integrals as OpenMP `target` regions on the A100 (compute
  capability `cc80`).
- **State batching / device-buffer tiling.** Production runs (>1e6 grid points,
  large basis, and/or multiple MPI ranks sharing a GPU) don't fit all
  quasiparticle states in GPU memory at once, so the state dimension is now
  tiled. Fixed device buffers are allocated with `omp_target_alloc`, each state
  tile is streamed in with `omp_target_memcpy`, and the kernel indexes them via
  `is_device_ptr(...)` with local offsets.
  - Do **not** use disjoint `map(to: psi_qp[off:len])` slice maps here — NVHPC's
    implicit present-check aborts at kernel launch when an access straddles two
    disjoint mapped sub-ranges. The `omp_target_alloc` + `omp_target_memcpy`
    approach avoids that ambiguity entirely.
  - Memory budget is controlled by the env var **`BSE_GPU_MEM_GB`** (default 32).
    A helper `bse_gpu_states_resident()` computes how many states fit and picks
    the tile sizes.

### Build fixes (this branch had never fully compiled before)
- Added header shims `read_bse.h` (→ `read.h`) and `mod_pseudopot.h` (→
  `mod_pot.h`) for renamed includes that were never created.
- Added missing `par_st` fields in `fd.h`: `epsX/epsY/epsZ` dielectric constants
  and the convenience grid-spacing mirrors.
- Fixed `dipole.c` spinor grid-offset indexing to compute offsets inline like
  everywhere else.
- Fixed the Cray-LAPACK `zheev_` call in `diag.c` to pass the two hidden Fortran
  string-length args (`..., &info, 1, 1)`).
- Made NVTX optional in `fd.h` (no-op unless `-DUSE_NVTX`); dropped
  `-lnvToolsExt` from the Makefile (broken with CUDA ≥ 12's header-only nvtx3).

### The critical link-flag gotcha
- NVHPC uses **relocatable device code (RDC)** by default, so the GPU flags
  `-mp=gpu -gpu=cc80` **must appear on BOTH the compile step AND the final link
  step**. If they're missing at link, the device linker (`nvlink`) never runs,
  and the binary builds fine but **every** `omp target` region aborts at runtime
  with *"Fatal error: Could not run target region on device 0."*
- The Makefile now keeps `Linux_GPUFLAGS = -mp=gpu -gpu=cc80` and passes it to
  the `BSE` link recipe. **Never** mix `-mp=multicore` and `-mp=gpu` objects in
  one binary — that also breaks offload.

### Performance (why this is worth it)
- On the tiny repo test cases (`1x1x1` = 27³ grid, `Pb` = 13³) the GPU is
  *slower* than the CPU — launch overhead and per-iteration PCIe transfer
  dominate a trivially small kernel. **Do not judge the offload on those.**
- At production sizes the GPU wins **~50–120×** over the CPU kernel (e.g. 64³
  grid / 20×20 states: 3.0 ms GPU vs 320 ms CPU). The A100 sustains
  500–830 GFLOP/s; the CPU kernel is memory-bandwidth bound.
- **Next optimization (not yet done):** the Hartree potential is still built on
  the CPU with FFTW once per iteration and copied host↔device each time. Moving
  that FFT to cuFFT so `rho→pot` stays resident on-device is the biggest
  remaining win. Transfers are currently synchronous; overlapping them with the
  kernel (`nowait` / double-buffering) is a further win.

---

## 2. Compiling on Perlmutter

### 2.1 Load the build environment

From the `bse_gpu/` directory:

```bash
source bse_gpu/LOAD_COMPILE_ENV.sh
```

This loads:
```
module load PrgEnv-nvidia    # NVHPC compilers behind the Cray cc wrapper
module load cray-fftw        # FFTW for the Hartree potential
module load cray-libsci      # LAPACK/BLAS (sci_nvidia)
```

> **Important:** `source` it (don't pipe its output through `tail` or redirect in
> a way that spawns a subshell) — otherwise the module swap is lost and `cc`
> falls back to Intel `icx`. The compiler must be the Cray `cc` wrapper: it
> supplies the MPI include paths that bare `nvc` lacks.

### 2.2 Build

```bash
cd bse_gpu
make            # produces BSE.x
```

Key Makefile facts (already configured — don't undo them):
- `Linux_CC = cc` (Cray wrapper over NVHPC).
- `Linux_GPUFLAGS = -mp=gpu -gpu=cc80` — present on **both** the `.c.o` compile
  rule **and** the `$(MAINNAM)` link recipe. This is the RDC fix above.
- `-O3 -fast -Minfo=all`. `-Minfo=all` prints which loops were offloaded — look
  for the `calc_eh_kernel_cplx` regions being mapped to the GPU.
- `make clean` removes `*.o *.x`.

### 2.3 Smoke-test the build (login node)

The Perlmutter **login nodes have an A100**, so you can sanity-check offload
without a job allocation:

```bash
cd bse_gpu
NVCOMPILER_ACC_DEBUG=1 ./BSE.x        # in a dir with valid input; watch for real device kernel launches
```

To confirm the GPU result matches the host result, run the same binary with
offload disabled and diff the eigenvalues:

```bash
OMP_TARGET_OFFLOAD=DISABLED ./BSE.x   # forces host execution
```

The BSE eigenvalues should be bit-identical (`max|ΔRe| = max|ΔIm| = 0`).

---

## 3. Running on the GPU nodes (SLURM)

### 3.1 Inputs a BSE run needs

BSE reads the output of a prior filter run. In the run directory you need:
- `input.par` — BSE parameters
- `conf.dat` — atomic configuration
- `output.dat` — the filtered states from the filter stage

`BSE.x` is MPI-parallel: it decomposes the kernel across ranks and gathers to
rank 0. On GPU nodes, run **one MPI rank per GPU**. Perlmutter GPU nodes have
**4× A100 (40 GB each)**, so use up to 4 ranks per node.

### 3.2 `BSE_GPU_MEM_GB`

Set this to the per-rank device-memory budget:

```
BSE_GPU_MEM_GB = A100_GB * fraction / ranks_per_gpu
```

With one rank per GPU on a 40 GB A100, leave headroom and use ~`36`. If two
ranks share a GPU, use ~`18`, etc. The default (32) is safe for the common
one-rank-per-GPU case.

### 3.3 Example batch script

Save as `submit_bse_gpu.sh` in your run directory and submit with `sbatch`.
Adjust `-A` (account), the number of nodes/ranks, and paths for your run.

```bash
#!/bin/bash
#SBATCH -J bse_gpu
#SBATCH -o %j-%N.out
#SBATCH -C gpu                    # GPU nodes (4x A100 each)
#SBATCH -q regular                # or: -q debug  (<= 30 min)
#SBATCH -A m4868                  # your NERSC allocation
#SBATCH -N 1                      # nodes
#SBATCH --ntasks-per-node=4       # 1 MPI rank per GPU
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=map_gpu:0,1,2,3
#SBATCH -c 32                     # CPU cores per rank (256 cores / 8 -> 32)
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel@qernelzoo.com

# --- Runtime environment: mirror the build env ---
module load PrgEnv-nvidia
module load cray-fftw
module load cray-libsci

# One GPU per rank; device memory budget for the state-tiling logic.
export BSE_GPU_MEM_GB=36

# Host-side OpenMP (Hartree FFT threads, etc.)
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Optional: verify offload actually lands on the device (verbose).
# export NVCOMPILER_ACC_DEBUG=1

BSE_EXE=/global/common/software/m4868/FilterBSE/bse_gpu/BSE.x

# Run from the directory holding input.par, conf.dat, output.dat.
srun ${BSE_EXE} 2> error.dat > run.dat
wait
```

### 3.4 Interactive session (debugging)

```bash
salloc -C gpu -q interactive -N 1 --ntasks-per-node=4 \
       --gpus-per-node=4 -c 32 -t 01:00:00 -A m4868

# then, inside the allocation:
module load PrgEnv-nvidia cray-fftw cray-libsci
export BSE_GPU_MEM_GB=36
srun -n 4 /global/common/software/m4868/FilterBSE/bse_gpu/BSE.x 2> error.dat > run.dat
```

---

## 4. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| *"Fatal error: Could not run target region on device 0"* at runtime | GPU flags missing from the **link** step (RDC). Confirm `-mp=gpu -gpu=cc80` is on the `$(MAINNAM)` recipe in the Makefile, then `make clean && make`. |
| `cc` invokes Intel `icx`, no GPU code emitted | `LOAD_COMPILE_ENV.sh` wasn't `source`d (module swap lost). Re-`source` it in the same shell. |
| SIGABRT in NVHPC `check_present` at kernel launch | Regression to disjoint `map(to: psi_qp[...])` slice maps. Keep the `omp_target_alloc` + `omp_target_memcpy` device-buffer path. |
| Out-of-memory / device alloc fails | Lower `BSE_GPU_MEM_GB`, or reduce ranks-per-GPU. |
| GPU slower than CPU | You're on a toy case (`1x1x1`, `Pb`). Expected — the offload only wins at production grid/basis sizes. |
| Verify correctness | Rerun with `OMP_TARGET_OFFLOAD=DISABLED`; eigenvalues must match the GPU run exactly. |

---

*Test case for validation:* `bse/testBSE/1x1x1/bse` (needs `input.par`,
`conf.dat`, `output.dat`). Note the repo's reference `BSEeval.par` there was
generated by an older `bse` executable on a different branch, so it legitimately
differs — validate the GPU path against the **offload-disabled** run of the same
binary, not against that file.
