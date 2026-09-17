# Job options

`bse_cplx` ships with **two builds** that produce two different executables. All
configuration is via `input.par` (and one environment variable, `BSE_GPU_MEM_GB`);
the only command-line argument is `initUnsafe`.

## Building

| Makefile | Executable | Toolchain | Capabilities |
|---|---|---|---|
| `Makefile` (default) | `bse_cplx_gpu.x` | NVHPC + cray-libsci + cray-fftw | GPU offload (`USE_GPU_OFFLOAD`), distributed ScaLAPACK (`USE_SCALAPACK`), libsci (`USE_LIBSCI`) |
| `Makefile_mkl` | `bse_cplx.x` | Intel + MKL + FFTW | CPU-only, single-node dense solve (no GPU, no ScaLAPACK) |

=== "GPU + distributed (default)"

    ```sh
    source LOAD_COMPILE_ENV_NVHPC.sh    # module load PrgEnv-nvidia cray-fftw cray-libsci
    make                                # -> bse_cplx_gpu.x
    ```

    Key flags: `-DFFTW_FFT -DUSE_LIBSCI -DUSE_GPU_OFFLOAD -DUSE_SCALAPACK`, with
    `-mp=gpu -gpu=cc80` on **both** compile and link (RDC requirement), `-O3
    -fast`. cray-libsci (the threaded `libsci_nvidia_mpi_mp`, which supplies
    `pzheevd`) is auto-linked by the `cc` wrapper — it is deliberately **not**
    linked with an explicit `-lsci` to avoid a double-link warning.

=== "CPU-only (MKL)"

    ```sh
    source LOAD_COMPILE_ENV.sh          # Intel environment
    make -f Makefile_mkl               # -> bse_cplx.x
    ```

    This build has no GPU offload and no ScaLAPACK, so it runs the serial dense
    `LAPACKE_zheev` eigensolve on a single node.

!!! warning "Switching toolchains needs a clean"
    The two builds share `.o` object names, so run `make clean` when switching
    between them. Loading modules: source the env script **in-shell** (no pipes —
    a pipe spawns a subshell and the `PrgEnv` swap is silently lost, leaving `cc`
    pointed at Intel).

Other historical Makefiles exist (`Makefile_nv`, `Makefile_gpu`, `Makefile.save`)
but the documented, current pair is `Makefile` and `Makefile_mkl`.

## Launching

`bse_cplx` uses **one MPI rank per GPU**, with the direct and exchange kernels
split across ranks by parity (see [Parallelization](parallelization.md)). A
representative launch:

```sh
export OMP_NUM_THREADS=32
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export BSE_GPU_MEM_GB=36

srun -n 4 --ntasks-per-node=4 --gpus-per-node=4 \
     --gpu-bind=map_gpu:0,1,2,3 -c 32 \
     ./bse_cplx_gpu.x > run.dat 2> error.dat
```

!!! tip "Use an even number of ranks ≥ 2"
    The Coulomb kernel splits work so that **even ranks compute the direct term
    and odd ranks compute the exchange term concurrently**. A single rank does
    not exercise this split correctly — run with ≥ 2 ranks on the distributed
    build, or use the MKL serial build for a one-node CPU run. See
    [Parallelization](parallelization.md).

The run directory must contain `input.par` and the filter output (`output.dat`
for the safe path, or `unsafe_input.par` + `eval.par` + `conf.par` + `psi.par`
for the unsafe path).

## Run modes

### Safe vs unsafe input

Set by the first CLI argument (`argv[1]`): absent/`0` reads `output.dat`; `1`
reads the `.par` files. See [Input files](input-files.md#two-input-paths).

### GPU on/off

`gpuAccel` in `input.par` (default `1`). At runtime the GPU path is used only if
the binary was built with `-DUSE_GPU_OFFLOAD` **and** `omp_get_num_devices() > 0`.
Setting `OMP_TARGET_OFFLOAD=DISABLED` transparently forces the host path — handy
for an A/B correctness check (the GPU and CPU kernels are bit-identical):

```sh
# run once each, then diff BSEeval.par
OMP_TARGET_OFFLOAD=DISABLED srun ... ./bse_cplx_gpu.x
srun ... ./bse_cplx_gpu.x
```

### Distributed vs serial

Automatic. On the ScaLAPACK build, `mpi_size > 1` selects the block-cyclic kernel
+ distributed `pzheevd`; `mpi_size == 1` (or the MKL build) uses the full-matrix
kernel + serial `LAPACKE_zheev`.

### Partial / staged runs

| Flag | Effect |
|---|---|
| `noCalcExciton = 1` | stop after the single-particle (QP) dipole/spin properties |
| `calcCoulombOnly = 1` | exit after building the e–h kernel (before the BSE solve) |
| `coulombDone = 1` | skip the kernel build and load `direct.dat` / `exchange.dat` from disk |

These let you split a large job: compute the expensive kernel once, then re-run
the (cheaper) diagonalization and analysis without redoing it.

## GPU memory budget: `BSE_GPU_MEM_GB`

The one environment variable the code reads. It is the **per-rank / per-GPU device
memory budget in GB** for the resident quasiparticle state set (default `36.0`,
with ~10% headroom reserved). It drives automatic state-tiling in both the direct
and exchange kernels:

- If the required states fit the budget, the kernel takes the fast single-map
  path.
- If not, it **automatically tiles** the states (or falls back to CPU for that
  kernel). There is no correctness impact — only speed.

Set it to `A100_GB × fraction / ranks_per_gpu`. On a 40 GB A100 with one rank per
GPU, `36` is the recommended value.

!!! note "Tiling has a Hartree-recompute cost"
    The tiled exchange kernel recomputes the host Hartree FFT once per
    (hole-tile × electron-tile) pass, so pass count matters. The code picks tile
    shapes that minimize the number of passes, but the takeaway is: give each GPU
    as much budget as fits so tiling is avoided when possible.

## Restarting and checkpointing

The Coulomb kernel is the expensive part, so restart support centers on it:

| Setting | Reads | Behavior |
|---|---|---|
| `coulombDone = 1` | `direct.dat` / `exchange.dat` | load the whole kernel, skip straight to the BSE solve |
| `restartCoulomb = 1` | `direct-<rank>.dat` / `exchange-<rank>.dat` | resume a partially-computed kernel from per-rank dumps |
| `saveCheckpoints = 1` | — | write job checkpoints |
| `restartFromChk = <id>` | checkpoint `<id>` | restart from a checkpoint |

!!! info "Distribute-on-read"
    On the distributed build, `coulombDone = 1` loads only each rank's own
    block-cyclic tile of the kernel (not the full matrix), so restart preserves
    the memory savings.
