# bse_cplx with optional GPU offload

`bse_cplx` keeps its **native `double complex`** typing everywhere. The Coulomb
kernel (`coulomb.c : calc_eh_kernel_cplx` — the direct `K^d` and exchange `K^x`
electron–hole integrals) can now be offloaded to an NVIDIA GPU via OpenMP
`target teams distribute`. Everything else (filter I/O, Hartree FFT, the dense
`zheev` diagonalization) still runs on the host. The GPU path is **selected at
runtime** and produces results **bit-identical** to the CPU path.

## Two builds (they coexist)

| Makefile | Toolchain | Output | Notes |
|---|---|---|---|
| `Makefile` (default) | NVHPC + cray-libsci, `-mp=gpu -gpu=cc80` | `bse_cplx_gpu.x` | GPU offload compiled in (`-DUSE_GPU_OFFLOAD`) |
| `Makefile_mkl` | Intel + MKL (original, verbatim) | `bse_cplx.x` | CPU only; the paradigm `bse_real` also uses |

The two share `.o` names, so run **`make clean`** when switching toolchains.

### Build the GPU version

```bash
# Load the env in THIS shell -- do NOT pipe module output anywhere: a pipe
# spawns a subshell and the swap is lost, leaving cc as Intel icx.
module load PrgEnv-nvidia
module load cray-fftw cray-libsci
make clean && make            # -> bse_cplx_gpu.x
```

`-mp=gpu -gpu=cc80` is on **both** the compile and the link step. NVHPC uses
relocatable device code (RDC), so omitting the flags at link builds a binary
whose every `omp target` region aborts at runtime with *"Could not run target
region on device 0"*.

### Build the original CPU version

```bash
module load PrgEnv-intel cray-fftw   # the default environment
make clean && make -f Makefile_mkl bse_cplx   # -> bse_cplx.x
```

## Toggling the GPU at runtime

Add to `input.par` (default is `1`):

```
gpuAccel = 1     # 1 = use GPU when a device is available; 0 = force CPU
```

The choice is gated by `omp_get_num_devices()`, so it degrades gracefully:

- No GPU on the node, or `gpuAccel = 0` → runs on the CPU.
- `OMP_TARGET_OFFLOAD=DISABLED` → forces the host path (handy for A/B checks).
- `BSE_GPU_MEM_GB` (default 32) budgets the resident quasiparticle-state set;
  if a kernel's states would exceed it, that kernel falls back to the CPU.

The run log prints which path was taken, e.g.
`Coulomb kernel: GPU offload ENABLED (1 device(s) visible)`.

## Verifying correctness

Run the same binary twice and diff the eigenvalues — they must match exactly:

```bash
gpuAccel=1 ...                          # GPU
OMP_TARGET_OFFLOAD=DISABLED ...         # same binary, host path
diff run_gpu/BSEeval.par run_cpu/BSEeval.par   # -> identical
```

Validated on a login-node A100 (c_2x2x2 filter output, 4h/4e): `BSEeval.par`,
`bsRE.dat`, `bsIM.dat` all bit-identical between GPU and CPU, `max|Δ| = 0`.

## How the offload is structured (native complex, no split-double)

- Hole states (direct) or hole+electron states (exchange) are mapped **once** to
  the device with a single contiguous `target enter data map` — no
  `omp_target_alloc` device pointers, and no disjoint slice maps (which trip
  NVHPC's present-check). This is the lesson from `bse_gpu` applied cleanly.
- Each `(a,b)`/`(a,i)` iteration recomputes the Hartree potential on the host
  (FFTW) and pushes only `pot_htree` to the device; one GPU team handles each
  hole/e–h pair and reduces the grid integral across its threads. Only a small
  result block returns to the host, so the full `direct`/`exchange` matrices
  never leave host memory.
- Native `double complex` arithmetic (`conj`, complex `*`) runs on the device
  (NVHPC ≥ 26). OpenMP C has no complex reduction, so the per-grid sum reduces
  on split real/imag doubles and is recombined — the only concession to the GPU.

## On the MKL → libsci switch

The only MKL use was the dense `zheev` in `diag.c`. Under cray-libsci we keep
`LAPACKE_zheev(LAPACK_ROW_MAJOR, …)` verbatim (libsci ships `lapacke.h`), so
there is **no eigenvector-layout change** and no downstream impact. `fd.h`
selects `<lapacke.h>` vs `<mkl.h>` on `USE_LIBSCI`, and `diag.c` uses
`omp_set_num_threads` instead of `mkl_set_num_threads` there. The eigensolve
stays on the host — appropriate for the typical ≤ 40000×40000 BSE matrix.
