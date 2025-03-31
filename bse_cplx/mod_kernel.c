#include "mod_kernel.h"

void mod_kernel(
  double complex*   psi_qp,
  double complex**  direct,
  double complex**  exchange,
  double complex*   pot_bare,
  double complex*   pot_screened,
  index_st*         ist,
  par_st*           par,
  flag_st*          flag,
  parallel_st*      parallel
  ){
  
  const int         mpir = parallel->mpi_rank;
  
  if (mpir == 0){
    write_separation(stdout, "T");
    printf("\n3.\tCOMPUTING ELEC-HOLE INTERACTION KERNEL | %s\n", get_time());
    write_separation(stdout, "B"); fflush(stdout);
  } 
  
  if (mpir == 0) printf("\nThe number of electron-hole pairs in the exciton basis = %ld\n", ist->n_xton);

  ALLOCATE(direct,   ist->n_xton * ist->n_xton, "direct");
  ALLOCATE(exchange, ist->n_xton * ist->n_xton, "exchange");
    
  if (0 == flag->coulombDone){
    if (mpir == 0) printf("Computing complex e-h kernel\n"); 
    fflush(0);

    calc_eh_kernel_cplx(
      psi_qp, pot_bare, pot_screened, *direct, *exchange,
      ist, par, flag, parallel
    );
    
  } else if (1 == flag->coulombDone){
    int done_flag;
    long a, i, b, j;
    printf("Loading in Coulomb matrix elements from files\n"); 
    fflush(0);

    done_flag = load_coulomb_mat(*direct, "direct.dat", &a, &b, &i, &j, ist);
    if (done_flag != 1){
      printf("ERROR: direct matrix elements not done!\n");
      exit(EXIT_FAILURE);
    }

    done_flag = load_coulomb_mat(*exchange, "exchange.dat", &a, &b, &i, &j, ist);
    if (done_flag != 1){
      printf("ERROR: exchange matrix elements not done!\n");
      exit(EXIT_FAILURE);
    }
  }
  
  if (flag->calcCoulombOnly == 1){
    printf("Exitng program after computing Coulomb matrix elements | %s\n", get_time());
    exit(0);
  }

  return;
}