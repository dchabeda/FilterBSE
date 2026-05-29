#include "fd.h"
#include "Hmat.h"
#include "energy.h"
#include "ortho.h"
#include "filter.h"
#include "aux.h"

// Periodic path: on each k-group master, gather -> orthogonalize -> diagonalize ->
// sigma_E -> write eval-[k].dat for every k-point owned by this group. Replaces the
// separate gather/ortho/diag/sigma/output modules when flag->periodic is set.
void run_periodic_postfilter(
    double *psi_rank,
    double *pot_local,
    xyz_st *R,
    grid_st *grid,
    vector *G_vecs,
    vector *k_vecs,
    zomplex *LS,
    nlc_st *nlc,
    long *nl,
    double *ksqr,
    double *eig_vals,
    double *sigma_E,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);

void mod_diag(
    double *psitot,
    double *pot_local,
    double *eig_vals,
    double *sigma_E,
    grid_st *grid,
    vector *G_vecs,
    vector *k_vecs,
    zomplex *LS,
    nlc_st *nlc,
    long *nl,
    zomplex *an,
    double *zn,
    double *ene_targets,
    double *ksqr,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel);
