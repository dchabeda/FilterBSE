# Periodic filter-diagonalization example: 1x1x1 CsPbI3
#
# This exercises the k-point-grouped periodic path of Filter_periodic.x
# (build it in ../../  ->  ../../../periodic_filter via `make`).
#
# Inputs
# ------
# input.par          : standard filter input, plus the periodic switches:
#                        periodic = 1            turn on the periodic path
#                        box_z = 11.884          length (Bohr) of the periodic
#                                                z box (= lattice c). x,y auto-size
#                                                to the cluster (1D/2D periodicity).
#                      The periodic kinetic operator kinetic_k now caps 0.5|k+G|^2 at
#                      KEmax (default 20), matching the non-periodic path, so the
#                      spectrum is bounded and the default iterative energy range
#                      (approxEnergyRange = 0) is stable. No KEmax/approxEnergyRange
#                      tuning is needed for this example.
# periodic_input.par : lattice vectors a1/a2/a3, lengths a/b/c, k-grid nk1/nk2/nk3,
#                        readKPath (0 = uniform grid, 1 = read kpath.par),
#                        nBands.
#
# k-point MPI grouping
# --------------------
# The world communicator is split into k-groups. With NK = nk1*nk2*nk3 k-points
# and NF = nFilterCycles:
#   mpi_size >= NK :  group_size = mpi_size/NK ranks per k (mpi_size must be a
#                     multiple of NK; group_size must divide NF). Each group does
#                     one k-point; filtering of NF cycles is split across the group.
#   mpi_size <  NK :  one group = all ranks, looping over all k-points; NF is split
#                     across all ranks (mpi_size must divide NF).
# Each k-group's master writes eval-<k>.dat (and psi-<k>.dat). The k-point equal to
# diag_k_idx (default 0) also writes output.dat for downstream BSE.
#
# Run (on a compute node)
# -----------------------
#   srun -n <NK or a multiple> -c <cpus> ./Filter_periodic.x | tee run.dat
# e.g. nk3 = 2  ->  srun -n 2 -c 16 ./Filter_periodic.x
#
# A quick single-rank smoke test (loops both k-points on one rank):
#   OMP_NUM_THREADS=8 ./Filter_periodic.x | tee run.dat
#
# Output: eval-0.dat, eval-1.dat, ...  (columns: index  eigenvalue[a.u.]  sigma_E)
# Concatenate eval-<k>.dat across k in order to assemble a band structure.
