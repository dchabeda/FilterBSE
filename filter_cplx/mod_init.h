#include "fd.h"
#include "read.h"
#include "init.h"
#include "init_periodic.h"
#include "aux.h"

void mod_init(
    grid_st*      grid,
    double**      gridx,
    double**      gridy,
    double**      gridz,
    xyz_st**      R,
    atom_info**   atom,
    double**      ene_targets,
    double**      ksqr,
    lattice_st**  lattice,
    vector**      G_vecs,
    vector**      k_vecs,
    index_st*     ist,
    par_st*       par,
    flag_st*      flag,
    parallel_st*  parallel
);



