#include "fd.h"

void read_filter_input(
    flag_st*      flag,
    grid_st*      grid,
    index_st*     ist,
    par_st*       par,
    parallel_st*  parallel);

void read_conf(
    xyz_st*       R,
    atom_info*    atm, 
    index_st*     ist, 
    par_st*       par, 
    flag_st*      flag, 
    parallel_st*  parallel
);

void read_pot(pot_st *pot, xyz_st *R, atom_info *atom, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void read_pot_file(FILE *pf, pot_st *pot, long j, long n, char *req);
void interpolate_pot(xyz_st *R, atom_info *atom, index_st *ist, par_st *par, parallel_st *parallel);
// void read_periodic_input(lattice_st *lattice, index_st *ist, par_st *par, flag_st *flag, parallel_st *parallel);
void calc_geom_par(xyz_st *R,atom_info *atm, index_st *ist, parallel_st *parallel);
double calc_bond_angle(long index1,long index2,long index3, xyz_st *R, parallel_st *parallel);
long assign_atom_number(char atyp[4]);
void assign_atom_type(char *atype,long j);
long get_number_of_atom_types(atom_info *atm,index_st *ist, long *list);
int assign_crystal_structure(char *crystal_structure);
int assign_outmost_material(char *outmost_material);
double get_ideal_bond_len(long natyp_1, long natyp_2, int crystalStructureInt);
  