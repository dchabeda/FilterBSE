#include "fd.h"

void main(int argc, char *argv[]){
	if (argc!=2) {printf("Usage: convert ngrid\n"); exit(EXIT_FAILURE);}
	long ngrid = atol(argv[1]);
	long nspinngrid = 2*ngrid;

	FILE *ppsi, *peval;
	double *psi, *eval, *deval;
	zomplex *zpsi;
	long i, j, k;

	//open and count number of lines in eval
	peval = fopen("eval.par" , "r");
	if (peval ==NULL) {
		printf("no eval in cwd\n");
		exit(EXIT_FAILURE);
	}
	for (i = j = 0; j != EOF; i++){
		j = fscanf(peval, "%*d %*g %*g");
	}
	fclose(peval);
	printf("Read %d states from eval\n",i);

	//allocate memory for eval and psi
	eval = (double*)  calloc(i,sizeof(double));
	deval= (double*)  calloc(i,sizeof(double));
	psi  = (double*)  calloc(ngrid,sizeof(double));
	zpsi = (zomplex*) calloc(2*i*nspinngrid,sizeof(zomplex));
	if(!zpsi){printf("insuffeciet memeory to convert to spinor!\n"); exit(EXIT_FAILURE);}

	ppsi = fopen("psi.par" , "r");
    if (ppsi ==NULL) {
        printf("no psi.par in cwd\n");
        exit(EXIT_FAILURE);
    }
    peval = fopen("eval.par" , "r");
    struct timeval tv;
    struct timezone tpz;
    gettimeofday(&tv,&tpz);
    srandom(tv.tv_sec);
    for(j = 0; j<i; j++){
    	fscanf(peval, "%*d %lg %lg", &eval[j], &deval[j]);
    	fread(psi, sizeof(double), ngrid, ppsi);
    	double coeff  = (double)rand() / (double)RAND_MAX;
    	for(k=0;k<ngrid;k++){
    		//spin + version
    		zpsi[2*j*nspinngrid+k].re = sqrt(coeff)*psi[k]; 
    		zpsi[2*j*nspinngrid+ngrid+k].re = sqrt(1-coeff)*psi[k];
    		//spin - version
    		zpsi[(2*j+1)*nspinngrid+k].im = sqrt(1-coeff)*psi[k];
    		zpsi[(2*j+1)*nspinngrid+ngrid+k].im = -sqrt(coeff)*psi[k];  
    	}
    }
    fclose(ppsi);

    ppsi = fopen("psi.out" , "w");
    fwrite (zpsi,sizeof(zomplex),2*i*nspinngrid,ppsi);
    fclose(ppsi);
    peval= fopen("eval.out", "w");
	for (j = 0; j < i; j++) {fprintf (ppsi,"%ld %.16g %g\n", j, eval[j], deval[j]); fprintf (ppsi,"%ld %.16g %g\n", j, eval[j], deval[j]);}
  	fclose(peval);  


}