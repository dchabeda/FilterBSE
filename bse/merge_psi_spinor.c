#include "fd.h"

int countlines(const char *filename);

int main(int argc, char *argv[]){
	if (argc!=9) {printf("Usage: merge_psi_spinor.x homoeval homopsi lumoeval lumopsi ngrid nhomo nlumo deps\n"); exit(EXIT_FAILURE);}
	long ngrid = atol(argv[5]);
	long nspinngrid = 2*ngrid;
	FILE *pevalhomo,  *pevallumo, *ppsihomo, *ppsilumo, *peval, *ppsi;
	zomplex *psihomo, *psilumo;
	double *eval, *deval;
	double evalloc, deloc;
	double fermiEnergy = -0.19;
	long i, a, ieof, j, k, ihomo, ilumo;
	long nhomo = atol(argv[6]);
	long nlumo = atol(argv[7]);
	double deps = atof(argv[8]);
	int nlines_ho, nlines_lu;
	
	// Count the number of lines in the eval.dat files
	nlines_ho = countlines(argv[1]);
	nlines_lu = countlines(argv[3]);

	printf("\nRunning merge_psi.c with the following inputs:\n");
	printf("homo eval file: %s, length: %d lines\n", argv[1], nlines_ho);
	printf("homo psi file: %s\n", argv[2]);
	printf("lumo eval file: %s, length: %d lines\n", argv[3], nlines_lu);
	printf("lumo psi file: %s\n", argv[4]);
	printf("ngrid: %ld\n", ngrid);
	printf("nspinngrid: %ld\n", nspinngrid);
	printf("nhomo = %ld; nlumo = %ld\n", nhomo, nlumo);
	printf("deps = %g\n", deps);
	printf("nhomo+nlumo = %ld\n", nhomo+nlumo);

	

	// allocate memory for the energies and variances
	if ((eval = (double *) calloc(nhomo+nlumo, sizeof(double))) == NULL) {
		printf("eval error"); 
		exit(1);
	}
	if ((deval = (double *) calloc(nhomo+nlumo, sizeof(double))) == NULL) {
		printf("deval error!"); 
		exit(1);
	}
	
	// Open homoeval; determine the index of the HOMO, ihomo.
	pevalhomo = fopen(argv[1],"r");
	if (pevalhomo == NULL){
		printf("no eval provided for hole states!\n");
		exit(EXIT_FAILURE);
	}
	// // determine index
	for (i = ieof = 0; ieof != EOF; i++){
		if (i >= nlines_ho) break; // break the loop if we iterate past the number of lines in the file
        ieof = fscanf(pevalhomo, "%ld %lg %lg", &a, &evalloc, &deloc);
        if (deloc < deps && evalloc < fermiEnergy) ihomo = i; // Grabs the last state with small variance before exceeding the fermiEnergy 
		if (evalloc > fermiEnergy) break; // do not search through the electron states
  	}
  	fclose(pevalhomo);
	printf("\nThe index of the homo is: %ld\n", ihomo);
	
	// Open lumoeval, determine the index of the LUMO, ilumo.
	pevallumo = fopen(argv[3],"r");
	if (pevallumo == NULL){
		printf("no eval provided for electron states!\n");
		exit(EXIT_FAILURE);
	}
	// // determine index
	for (i = ieof = 0; i != EOF; i++){
		ieof = fscanf(pevallumo, "%ld %lg %lg", &a, &evalloc, &deloc);
		//printf("The current i: %ld\n", i);fflush(0);
		if (deloc < deps && evalloc > fermiEnergy){
			ilumo = i;
			break; // Grabs the first state above the fermiEnergy with small variance
		} 
	}
	fclose(pevallumo);
	printf("\nThe index of the lumo is: %ld\n", ilumo); //fflush(0);
	
	// Populate *eval and *deval with the nhomo states below ihomo
	pevalhomo = fopen(argv[1], "r");
	
	printf("\nGathering hole states...\n"); fflush(0);
	for (i = ieof = 0; i <= ihomo; i++){
		if (i <= ihomo - nhomo){
			fscanf(pevalhomo, "%ld %lg %lg", &a, &evalloc, &deloc); // throw away these values of eval and de
		}
		if (i > ihomo - nhomo){
			fscanf(pevalhomo, "%ld %lg %lg", &a, &evalloc, &deloc); // store these values that are within nhomo of valence band edge
			eval[ieof] = evalloc; 
			deval[ieof] = deloc;
			printf("%lg %lg\n", eval[ieof], deval[ieof]);
			ieof++;	
		}
	}
	
	printf("Hole energies successfully acquired\n"); // fflush(0);
	
	// Populate *eval and *deval with the nlumo states above ilumo
    printf("\nGathering electron states...\n"); fflush(0);
	pevallumo = fopen(argv[3], "r");
	for (i = 0; i < ilumo; i++){
		fscanf(pevallumo, "%ld %lg %lg", &a, &evalloc, &deloc);
		printf("Discard state %ld: %lg\n", i, evalloc);
	} // Just to move the file pointer to the bottom of the electron states
	
	for (i = 0; i < nlumo; i++){
		fscanf(pevallumo, "%ld %lg %lg", &a, &eval[i+nhomo], &deval[i+nhomo]); // store these values that are within nhomo of valence band edge
        printf("%lg %lg\n", eval[i+nhomo], deval[i+nhomo]);
		}
	printf("Electron energies successfully acquired\n");

	// Print the new eval file with the merged evals
	printf("\nPrinting the new eval.dat\n"); //fflush(0);
	peval = fopen("eval_new.dat", "w");
	for (i = 0; i < nhomo + nlumo; i++){
		fprintf(peval, "%ld %.16f %g\n", i, eval[i], deval[i]);
	}
	fclose(peval);

	//printf("\nMerge psi done\n"); fflush(0);
	
	// Populate the first nspinngrid*nhomo positions of *psi with psihomo[ihomo-nhomo->nhomo]
	// open lumoeval, determine the index of the lumo, ilumo
	// Populate *eval with nlumo states with de < deps
	// Populate the next nspinngrid*nlumo positions of *psi with psilumo[ilumo->ilumo+nlumo]


	//allocate memory for psi
	if ((psihomo  = (zomplex*)  calloc(nspinngrid*nhomo,sizeof(zomplex))) == NULL) {
		printf("psihomo error\n");
		exit(1);
	};
	if ((psilumo  = (zomplex*)  calloc(nspinngrid*nlumo,sizeof(zomplex))) == NULL) {
                printf("psilumo error\n");
                exit(1);
        };
	ppsihomo = fopen(argv[2] , "r"); // get psi_homo
    if (ppsihomo ==NULL) {
        printf("no homo psi.par in cwd\n");
        exit(EXIT_FAILURE);
    }
    printf("\nGathering hole wavefunctions...\n"); fflush(0);
	fseek(ppsihomo, nspinngrid*(ihomo - nhomo + 1)*sizeof(zomplex), SEEK_SET);
	fread(psihomo, 1, sizeof(zomplex) * nspinngrid*nhomo, ppsihomo);
	fclose(ppsihomo);
	
	printf("\nGathering electron wavefunctions...\n"); fflush(0);
	ppsilumo = fopen(argv[4], "r"); // get psi_lumo
	fseek(ppsilumo, nspinngrid*ilumo*sizeof(zomplex), SEEK_SET);
	fread(psilumo, 1, sizeof(zomplex) * nspinngrid*nlumo, ppsilumo);
	fclose(ppsilumo);

	ppsi = fopen("psi_new.dat", "w");
	
	fwrite(psihomo, 1, sizeof(zomplex) * nspinngrid*nhomo, ppsi);
	fwrite(psilumo, 1, sizeof(zomplex) * nspinngrid*nlumo, ppsi);
	fclose(ppsi);
	
	printf("\nDone merging wavefunctions...\n"); fflush(0);
	return 0;

}

int countlines(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Could not open file");
        return -1;
    }

    int count = 0;
    char ch;

    // Read through the file character by character
    while ((ch = fgetc(file)) != EOF) {
        if (ch == '\n') {
            count++;
        }
    }

    fclose(file);
    return count;
}
