# bse_cplx

`bse_cplx` is the complex-valued Bethe–Salpeter equation (BSE) solver. It consumes
the quasiparticle (QP) states and energies produced by
[`filter_mpi`](../filter_mpi/index.md) and computes **excitons** and **optical
properties** of a semiconductor nanostructure: exciton energies, binding
energies, oscillator strengths, magnetic and rotational strengths, and spin /
angular-momentum character.

It is written in C with **native `double complex`** typing throughout, and offers
MPI distribution (ScaLAPACK), OpenMP threading, and optional GPU offload of the
Coulomb kernel.

## What it computes

<div class="grid cards" markdown>

- :material-database-import: **1. QP basis**

    Read the filter output, apply the \(\sigma_E\) / Fermi-energy window to pick
    converged hole and electron states, and load their wavefunctions into a
    node-shared buffer.

- :material-atom: **2. e–h kernel**

    Build the electron–hole interaction: the **direct** (screened Coulomb) and
    **exchange** (bare Coulomb) terms, via FFT convolutions.

- :material-grid: **3. BSE eigensolve**

    Assemble \(H = h_0 - (K^d + K^x)\) and diagonalize it for the **full**
    exciton spectrum (all eigenvectors — a physics requirement).

- :material-lightbulb-on: **4. Optical & angular**

    Contract the exciton coefficients with single-particle dipoles and
    angular-momentum operators to get spectra and state character.

</div>

## How to read these docs

| Page | What's in it |
|---|---|
| [Input files](input-files.md) | `input.par` reference, safe (`output.dat`) vs unsafe (`.par`) input paths |
| [Job options](job-options.md) | The two builds, launching, run modes, restart, GPU controls |
| [Codebase & module layout](codebase.md) | Every source file, `main.c` control flow, headers |
| [Pipeline](pipeline.md) | QP basis → kernel → BSE diag → optical/angular, and data flow |
| [Parallelization](parallelization.md) | MPI parity split, ScaLAPACK, node-shared psi, GPU auto-tiling, tips |
| [Outputs](outputs.md) | Every file the code writes and its format |

## Quick start

```sh
# 1. Load the GPU/ScaLAPACK build environment (Perlmutter)
source LOAD_COMPILE_ENV_NVHPC.sh     # PrgEnv-nvidia + cray-fftw + cray-libsci

# 2. Build the GPU + distributed binary
cd bse_cplx && make                  # produces bse_cplx_gpu.x

# 3. Run in a directory holding the filter output (output.dat) + input.par
export BSE_GPU_MEM_GB=36
srun -n 4 --gpus-per-task=1 --gpu-bind=closest -c 32 ./bse_cplx_gpu.x > run.dat 2> error.dat
```

For a CPU-only single-node build, use `Makefile_mkl` to produce `bse_cplx.x` — see
[Job options](job-options.md#building).

!!! info "Energies are in atomic units"
    Unless a column is explicitly labeled eV, all energies are in Hartree
    (\(\texttt{AUTOEV} = 27.2114\)).
