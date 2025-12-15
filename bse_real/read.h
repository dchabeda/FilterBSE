#include "fd.h"
#include "aux.h"

/**************************************************/

void read_input(
    flag_st*          flag, 
    grid_st*          grid, 
    index_st*         ist, 
    par_st*           par, 
    parallel_st*      parallel
);

/*************************************************/

void read_field(
    const char*     field, 
    const char*     key, 
    void*           var, 
    VarType         type, 
    const char*     tmp,
    int*            fk
);

/*************************************************/

void read_unsafe_input(
    double complex** psitot,
    double** eig_vals,
    double** sigma_E,
    xyz_st** R, 
    grid_st *grid,
    double** gridx,
    double** gridy,
    double** gridz,
    index_st *ist,
    par_st *par,
    flag_st *flag,
    parallel_st *parallel
);
  
/*************************************************/

void get_fmo_idxs(
    double*          eig_vals,
    double*          sigma_E,
    double           fermiE,
    double           secut,
    unsigned long    n_elems,
    unsigned long*   homo_idx
);

