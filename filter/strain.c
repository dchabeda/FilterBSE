#include "fd.h"
#include "vector.h"

/****************************************************************************/
/* read all neighbors of each atom */
void readNearestNeighbors(long nAtoms, int crystalStructure, vector *atomNeighbors, double *tetrahedronVolRef, int outmostMaterial) {
	
	FILE *pf;
	long iAtom, at_natyp, *neighbors_natyp; 
	char at[4], n1[4], n2[4], n3[4], n4[4];
	char strerror[200];
	int tmpi, numSe, numS, numCd, numAs, numIn, numP, numGa;
	double tmpx, tmpy, tmpz;
	double *bondLengths; 
	if ((neighbors_natyp = (long*)calloc(4,sizeof(double)))==NULL)nerror("neighbors_natyp");
	if ((bondLengths = (double*)calloc(4,sizeof(double)))==NULL)nerror("bondLengths");
	
	if ( access("allNeighborBonds.par", F_OK) != -1 ) {
		
		pf = fopen("allNeighborBonds.par", "r");
		
		for (iAtom = 0; iAtom < nAtoms; iAtom++) {
			fscanf(pf, "%s ", at);
			if ((! strcmp(at, "P1")) || (! strcmp(at, "P2")) || (! strcmp(at, "P3")) || (! strcmp(at, "P4"))) {
	fscanf(pf, "%i %lg %lg %lg %s %i %lg %lg %lg",
				 &tmpi, &tmpx, &tmpy, &tmpz,
				 n1, &tmpi, &atomNeighbors[4*iAtom].x, &atomNeighbors[4*iAtom].y, &atomNeighbors[4*iAtom].z);
	tetrahedronVolRef[iAtom] = 0.;
			}

			else {
	fscanf(pf, "%i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg",
				 &tmpi, &tmpx, &tmpy, &tmpz,
				 n1, &tmpi, &atomNeighbors[4*iAtom].x, &atomNeighbors[4*iAtom].y, &atomNeighbors[4*iAtom].z,
				 n2, &tmpi, &atomNeighbors[4*iAtom+1].x, &atomNeighbors[4*iAtom+1].y, &atomNeighbors[4*iAtom+1].z,
				 n3, &tmpi, &atomNeighbors[4*iAtom+2].x, &atomNeighbors[4*iAtom+2].y, &atomNeighbors[4*iAtom+2].z,
				 n4, &tmpi, &atomNeighbors[4*iAtom+3].x, &atomNeighbors[4*iAtom+3].y, &atomNeighbors[4*iAtom+3].z);
	
	// adjust positions of passivation ligands
	if (! strcmp(n3, "P1")) {
		strcpy(n3, n2);
		atomNeighbors[4*iAtom+2] = retScaledVector(atomNeighbors[4*iAtom+2], 1.0/0.55);
	}
	if (! strcmp(n4, "P1")) {
		strcpy(n4, n2);
		atomNeighbors[4*iAtom+3] = retScaledVector(atomNeighbors[4*iAtom+3], 1.0/0.55);
	}
	if (! strcmp(n3, "P2") && ! strcmp(n4, "P2")) {
		atomNeighbors[4*iAtom+2] = retScaledVector(atomNeighbors[4*iAtom+2], 1./0.25);
		atomNeighbors[4*iAtom+3] = retScaledVector(atomNeighbors[4*iAtom+3], 1./0.25);
		// Below is just for Eran's old geometry
		// atomNeighbors[4*iAtom+2] = retScaledVector(atomNeighbors[4*iAtom+2], 1./0.30);
		// atomNeighbors[4*iAtom+3] = retScaledVector(atomNeighbors[4*iAtom+3], 1./0.30);
	}
	else if (! strcmp(n4, "P2")) {
		atomNeighbors[4*iAtom+3] = retScaledVector(atomNeighbors[4*iAtom+3], 1.0/0.30);
		// Below is just for Eran's old geometry
		// atomNeighbors[4*iAtom+3] = retScaledVector(atomNeighbors[4*iAtom+3], 1.0/0.40);
	}
	
	at_natyp = assign_atom_number(at);
	neighbors_natyp[0] = assign_atom_number(n1);
	neighbors_natyp[1] = assign_atom_number(n2);
	neighbors_natyp[2] = assign_atom_number(n3);
	neighbors_natyp[3] = assign_atom_number(n4);
	
	for (int iNeighbor = 0; iNeighbor < 4; iNeighbor++) {
		bondLengths[iNeighbor] = 0.;
		// If neighbor is a passivation ligand, replace it with the corresponding semiconductor atom
		if ((neighbors_natyp[iNeighbor]==8) || (neighbors_natyp[iNeighbor]==9) || (neighbors_natyp[iNeighbor]==10) || (neighbors_natyp[iNeighbor]==11)) { 
			if ((outmostMaterial==0) && (at_natyp==0)) neighbors_natyp[iNeighbor]=7; // CdS, Center-Cd, Replace with S. 
			else if ((outmostMaterial==0) && (at_natyp==7)) neighbors_natyp[iNeighbor]=0; // CdS, Center-S, Replace with Cd. 
			else if ((outmostMaterial==0) && (at_natyp!=0) && (at_natyp!=7)) {
				sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp); 
				nerror(strerror);
			} 
			else if ((outmostMaterial==1) && (at_natyp==0)) neighbors_natyp[iNeighbor]=1; // CdSe, Center-Cd, Replace with Se. 
			else if ((outmostMaterial==1) && (at_natyp==1)) neighbors_natyp[iNeighbor]=0; // CdSe, Center-Se, Replace with Cd. 
			else if ((outmostMaterial==1) && (at_natyp!=0) && (at_natyp!=1)) {
				sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
				nerror (strerror);
			}
			else if ((outmostMaterial==2) && (at_natyp==2)) neighbors_natyp[iNeighbor]=16; // InP, Center-In, Replace with P. 
			else if ((outmostMaterial==2) && (at_natyp==16)) neighbors_natyp[iNeighbor]=2; // InP, Center-P, Replace with In. 
			else if ((outmostMaterial==2) && (at_natyp!=2) && (at_natyp!=16)) {
				sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
				nerror (strerror);
			}
			else if ((outmostMaterial==3) && (at_natyp==2)) neighbors_natyp[iNeighbor]=3; // InAs, Center-In, Replace with As. 
			else if ((outmostMaterial==3) && (at_natyp==3)) neighbors_natyp[iNeighbor]=2; // InAs, Center-As, Replace with In. 
			else if ((outmostMaterial==3) && (at_natyp!=2) && (at_natyp!=3)) {
				sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
				nerror (strerror);
			}
			else if ((outmostMaterial==4) && ((at_natyp==2)||(at_natyp==15))) neighbors_natyp[iNeighbor]=16; // Outmost: alloyInGaP; Center: In or Ga. Replace with P. 
			else if ((outmostMaterial==4) && (at_natyp==16)) neighbors_natyp[iNeighbor]=2; // Outmost: alloyInGaP; Center: P. Replace with In. This is a completely random choice. 
			else if ((outmostMaterial==4) && (at_natyp!=2) && (at_natyp!=15) && (at_natyp!=16)) {
				sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
				nerror (strerror);
			}
			else if ((outmostMaterial==5) && ((at_natyp==2)||(at_natyp==15))) neighbors_natyp[iNeighbor]=3; // Outmost: alloyInGaAs; Center: In or Ga. Replace with As. 
			else if ((outmostMaterial==5) && (at_natyp==3)) neighbors_natyp[iNeighbor]=2; // Outmost: alloyInGaAs; Center: As. Replace with In. This is a completely random choice. 
			else if ((outmostMaterial==5) && (at_natyp!=2) && (at_natyp!=15)) {
				sprintf(strerror,"Outmost layer is input as %d, but the surface is not cation terminated.\n", outmostMaterial);
				nerror (strerror);
			}
			else if ((outmostMaterial==6) && (at_natyp==15)) neighbors_natyp[iNeighbor]=3; // Outmost: GaAs; Center: Ga. Replace with As. 
			else if ((outmostMaterial==6) && (at_natyp==3)) neighbors_natyp[iNeighbor]=15; // Outmost: GaAs; Center: As. Replace with Ga. 
			else if ((outmostMaterial==6) && (at_natyp!=15) && (at_natyp!=3)) {
				sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
				nerror (strerror);
			}
			else if ((outmostMaterial==7) && (at_natyp==6)) neighbors_natyp[iNeighbor]=1; // ZnSe, Center-Zn, Replace with Se. 
			else if ((outmostMaterial==7) && (at_natyp==1)) neighbors_natyp[iNeighbor]=6; // ZnSe, Center-Se, Replace with Zn. 
			else if ((outmostMaterial==7) && (at_natyp!=1) && (at_natyp!=6)) {
				sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
				nerror (strerror);
			}
			else if ((outmostMaterial==8) && (at_natyp==6)) neighbors_natyp[iNeighbor]=7; // ZnS, Center-Zn, Replace with S. 
			else if ((outmostMaterial==8) && (at_natyp==7)) neighbors_natyp[iNeighbor]=6; // ZnS, Center-S, Replace with Zn. 
			else if ((outmostMaterial==8) && (at_natyp!=6) && (at_natyp!=7)) {
				sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
				nerror (strerror);
			}
			else if ((outmostMaterial==9) && ((at_natyp==3)||(at_natyp==16))) neighbors_natyp[iNeighbor]=15; // Outmost: alloyGaAsP; Center: As or P. Replace with Ga. 
			else if ((outmostMaterial==9) && (at_natyp==15)) neighbors_natyp[iNeighbor]=3; // Outmost: alloyGaAsP; Center: Ga. Replace with As. This is a completely random choice. 
			else if ((outmostMaterial==9) && (at_natyp!=3) && (at_natyp!=15) && (at_natyp!=16)) {
				sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
				nerror (strerror);
			}
                        else if ((outmostMaterial==11) && (at_natyp==15)) neighbors_natyp[iNeighbor]=16; // Outmost: GaP; Center: Ga. Replace with P. 
			else if ((outmostMaterial==11) && (at_natyp==16)) neighbors_natyp[iNeighbor]=15; // Outmost: GaP; Center: P. Replace with Ga. 
			else if ((outmostMaterial==11) && (at_natyp!=15) && (at_natyp!=16)) {
				sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
				nerror (strerror);
			}

		}
		bondLengths[iNeighbor] = retIdealBondLength(at_natyp, neighbors_natyp[iNeighbor], crystalStructure); 
		// printf("at_natyp = %d, iNeighbor = %d, neighbors_natyp[iNeighbor] = %d, bondLengths[iNeighbor]=%g \n", at_natyp, iNeighbor, neighbors_natyp[iNeighbor], bondLengths[iNeighbor]); 
	}
	tetrahedronVolRef[iAtom] = retRegularTetrahedronVolume(bondLengths[0], bondLengths[1], bondLengths[2], bondLengths[3]);
	// printf("The reference tetrahedron volume is %.7f \n", tetrahedronVolRef[iAtom]);  
			}
		}
		fclose(pf);
	}
	else {
		printf("\n\nNo allNeighborBonds.par file in current working directory -- the program is exiting!!!\n\n");
		fflush(stdout);
		exit(EXIT_FAILURE);
	}
	return;
}

/****************************************************************************/
/* compute strain scaling factor for each atom
 * strainScale = 1. + a4 * (Omega/Omega0 - 1),
 * where Omega is volume of tetrahedron formed by atom's nearest neighbors */

void calculateStrainScale(long nAtoms, double *tetrahedronVolRef, atm_st *atm, vector *atomNeighbors, double *a4Params, double *a5Params, double *strainScale) {

	FILE *pf;
	long iAtom;
	vector v1, v2, v3, v4;
	char at[4];
	double tmp2, tmp3, *tetrahedronVol, strain;
	int tmp1;

	if ((tetrahedronVol = (double *) calloc(nAtoms, sizeof(double))) == NULL) nerror("tetrahedronVol");

	if ( access("strain.par", F_OK) != -1 ) {
	
		printf("Reading strain scale from strain.par...\n");    
		pf = fopen("strain.par", "r");
		for (iAtom = 0; iAtom < nAtoms; iAtom++) {
			fscanf(pf, "%i %s %lg %lg %lg %lg %lg", &tmp1, at, &tmp2, &tmp3, 
					&tetrahedronVol[iAtom], &tetrahedronVolRef[iAtom], &strainScale[iAtom]);
		}
		fclose(pf);

	}
	else {
	
		// TODO: parallelize
		for (iAtom = 0; iAtom < nAtoms; iAtom++) {

			if ((atm[iAtom].natyp == 8) || (atm[iAtom].natyp == 9)) { 
				strainScale[iAtom] = 1.;
			}
			else {

				v1 = retSubtractedVectors(atomNeighbors[4*iAtom], atomNeighbors[4*iAtom+3]);
				v2 = retSubtractedVectors(atomNeighbors[4*iAtom+1], atomNeighbors[4*iAtom+3]);
				v3 = retSubtractedVectors(atomNeighbors[4*iAtom+2], atomNeighbors[4*iAtom+3]);

				tetrahedronVol[iAtom] = fabs(retDotProduct(v1, retCrossProduct(v2, v3)))/6.;
				strain = (tetrahedronVol[iAtom]/tetrahedronVolRef[iAtom] - 1.);
				strainScale[iAtom] = (1. + a4Params[atm[iAtom].natyp]*strain + a5Params[atm[iAtom].natyp]*cube(strain));
			}

		}
	}

	pf = fopen("strain.dat","w");
	for (iAtom = 0; iAtom < nAtoms; iAtom++) {
		fprintf(pf, "%ld %s %.8f %.8f %.8f %.8f %.8f\n", iAtom, atm[iAtom].atyp, a4Params[atm[iAtom].natyp], a5Params[atm[iAtom].natyp], 
				tetrahedronVol[iAtom], tetrahedronVolRef[iAtom], strainScale[iAtom]);
	}
	fclose(pf);

	free(tetrahedronVol);

	return;
}

/****************************************************************************/
double retRegularTetrahedronVolume(double bondLength1, double bondLength2, double bondLength3, double bondLength4) {
	
	vector v1, v2, v3, v4; 
	v1.x=sqrt(8./9.);   v1.y=0.;   v1.z=-1./3.;
	v2.x=-sqrt(2./9.);   v2.y=sqrt(2./3.);   v2.z=-1./3.;
	v3.x=-sqrt(2./9.);   v3.y=-sqrt(2./3.);   v3.z=-1./3.;
	v4.x=0.;   v4.y=0.;   v4.z=1.;
	
	v1 = retScaledVector(v1, bondLength1);
	v2 = retScaledVector(v2, bondLength2);
	v3 = retScaledVector(v3, bondLength3);
	v4 = retScaledVector(v4, bondLength4);

	vector side1 = retSubtractedVectors(v1, v4); 
	vector side2 = retSubtractedVectors(v2, v4); 
	vector side3 = retSubtractedVectors(v3, v4); 

	double regularTetrahedronVolume = fabs(retDotProduct(side1, retCrossProduct(side2, side3)))/6.;
	return regularTetrahedronVolume; 
}

/****************************************************************************/
