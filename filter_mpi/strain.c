#include "fd.h"


/****************************************************************************/
/* Read all neighbors of each atom */
void read_nearest_neighbors(vector *atom_neighbors, double *tetrahedron_vol_ref, long natoms, int crystal_structure, int outmost_material) {
	
	FILE *pf;
	long jatom, ctr_atom_natyp, *neighbors_natyp; 
	char ctr_atom[3], n1[3], n2[3], n3[3], n4[3];
	char strerror[200];
	int tmpi; // numSe, numS, numCd, numAs, numIn, numP, numGa;
	double tmpx, tmpy, tmpz;
	double *bond_lengths; 

	if ((neighbors_natyp = (long*) calloc(4, sizeof(double)))==NULL){
		fprintf(stderr, "OUT OF MEMORY: neighbors_natyp\n"); exit(EXIT_FAILURE);
	}
	if ((bond_lengths = (double*) calloc(4, sizeof(double)))==NULL){
		fprintf(stderr, "OUT OF MEMORY: bond_lengths\n"); exit(EXIT_FAILURE);
	}
	
	if ( access("allNeighborBonds.par", F_OK) != -1 ) {
		
		pf = fopen("allNeighborBonds.par", "r");
		
		for (jatom = 0; jatom < natoms; jatom++) {
			fscanf(pf, "%s ", ctr_atom);
			if ((! strcmp(ctr_atom, "P1")) || (! strcmp(ctr_atom, "P2")) || (! strcmp(ctr_atom, "P3")) || (! strcmp(ctr_atom, "P4"))){
				fscanf(pf, "%i %lg %lg %lg %s %i %lg %lg %lg",
				&tmpi, &tmpx, &tmpy, &tmpz,
				n1, &tmpi, &atom_neighbors[4*jatom].x, &atom_neighbors[4*jatom].y, &atom_neighbors[4*jatom].z);
				tetrahedron_vol_ref[jatom] = 0.;
			} 
      else {
        fscanf(pf, "%i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg",
              &tmpi, &tmpx, &tmpy, &tmpz,
              n1, &tmpi, &atom_neighbors[4*jatom].x, &atom_neighbors[4*jatom].y, &atom_neighbors[4*jatom].z,
              n2, &tmpi, &atom_neighbors[4*jatom+1].x, &atom_neighbors[4*jatom+1].y, &atom_neighbors[4*jatom+1].z,
              n3, &tmpi, &atom_neighbors[4*jatom+2].x, &atom_neighbors[4*jatom+2].y, &atom_neighbors[4*jatom+2].z,
              n4, &tmpi, &atom_neighbors[4*jatom+3].x, &atom_neighbors[4*jatom+3].y, &atom_neighbors[4*jatom+3].z);
        
	      // adjust positions of passivation ligands
        if (! strcmp(n3, "P1")) {
          strcpy(n3, n2);
          atom_neighbors[4*jatom+2] = retScaledVector(atom_neighbors[4*jatom+2], 1.0/0.55);
        }
        if (! strcmp(n4, "P1")) {
          strcpy(n4, n2);
          atom_neighbors[4*jatom+3] = retScaledVector(atom_neighbors[4*jatom+3], 1.0/0.55);
        }
        if (! strcmp(n3, "P2") && ! strcmp(n4, "P2")) {
          atom_neighbors[4*jatom+2] = retScaledVector(atom_neighbors[4*jatom+2], 1./0.25);
          atom_neighbors[4*jatom+3] = retScaledVector(atom_neighbors[4*jatom+3], 1./0.25);
          // Below is just for Eran's old geometry
          // atom_neighbors[4*jatom+2] = retScaledVector(atom_neighbors[4*jatom+2], 1./0.30);
          // atom_neighbors[4*jatom+3] = retScaledVector(atom_neighbors[4*jatom+3], 1./0.30);
        }
        else if (! strcmp(n4, "P2")) {
          atom_neighbors[4*jatom+3] = retScaledVector(atom_neighbors[4*jatom+3], 1.0/0.30);
          // Below is just for Eran's old geometry
          // atom_neighbors[4*jatom+3] = retScaledVector(atom_neighbors[4*jatom+3], 1.0/0.40);
        }
	
        ctr_atom_natyp = assign_atom_number(ctr_atom);
        neighbors_natyp[0] = assign_atom_number(n1);
        neighbors_natyp[1] = assign_atom_number(n2);
        neighbors_natyp[2] = assign_atom_number(n3);
        neighbors_natyp[3] = assign_atom_number(n4);
        
        for (int jneighbor = 0; jneighbor < 4; jneighbor++) {
          bond_lengths[jneighbor] = 0.;
          // If neighbor is a passivation ligand, replace it with the corresponding semiconductor atom
          if ((2 == neighbors_natyp[jneighbor]) || (3 == neighbors_natyp[jneighbor]) || (5 == neighbors_natyp[jneighbor])){ 
            if ((0 == outmost_material) && (48 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 16; // CdS, Center-Cd, Replace with S. 
            else if ((0 == outmost_material) && (16 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 48; // CdS, Center-S, Replace with Cd. 
            else if ((0 == outmost_material) && (ctr_atom_natyp != 48) && (ctr_atom_natyp != 16)) {
              sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmost_material, ctr_atom_natyp); 
              nerror(strerror);
            } 
            else if ((1 == outmost_material) && (48 == ctr_atom_natyp)) neighbors_natyp[jneighbor]= 34; // CdSe, Center-Cd, Replace with Se. 
            else if ((1 == outmost_material) && (34 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 48; // CdSe, Center-Se, Replace with Cd. 
            else if ((1 == outmost_material) && (ctr_atom_natyp != 48) && (ctr_atom_natyp != 34)) {
              sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmost_material, ctr_atom_natyp);
              nerror (strerror);
            }
            else if ((2 == outmost_material) && (49 == ctr_atom_natyp)) neighbors_natyp[jneighbor]= 15; // InP, Center-In, Replace with P. 
            else if ((2 == outmost_material) && (15 == ctr_atom_natyp)) neighbors_natyp[jneighbor]= 49; // InP, Center-P, Replace with In. 
            else if ((2 == outmost_material) && (ctr_atom_natyp != 49) && (ctr_atom_natyp != 15)) {
              sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmost_material, ctr_atom_natyp);
              nerror (strerror);
            }
            else if ((3 == outmost_material) && (49 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 33; // InAs, Center-In, Replace with As. 
            else if ((3 == outmost_material) && (33 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 49; // InAs, Center-As, Replace with In. 
            else if ((3 == outmost_material) && (ctr_atom_natyp != 49) && (ctr_atom_natyp != 33)) {
              sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmost_material, ctr_atom_natyp);
              nerror (strerror);
            }
            else if ((4 == outmost_material) && ((49 == ctr_atom_natyp)||(31 == ctr_atom_natyp))) neighbors_natyp[jneighbor] = 15; // Outmost: alloyInGaP; Center: In or Ga. Replace with P. 
            else if ((4 == outmost_material) && (15 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 49; // Outmost: alloyInGaP; Center: P. Replace with In. This is a completely random choice. 
            else if ((4 == outmost_material) && (ctr_atom_natyp != 49) && (ctr_atom_natyp != 31 ) && (ctr_atom_natyp != 15)) {
              sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmost_material, ctr_atom_natyp);
              nerror (strerror);
            }
            else if ((5 == outmost_material) && ((49 == ctr_atom_natyp)||(31 == ctr_atom_natyp))) neighbors_natyp[jneighbor] = 33; // Outmost: alloyInGaAs; Center: In or Ga. Replace with As. 
            else if ((5 == outmost_material) && (33 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 49; // Outmost: alloyInGaAs; Center: As. Replace with In. This is a completely random choice. 
            else if ((5 == outmost_material) && (ctr_atom_natyp != 49) && (ctr_atom_natyp != 31)) {
              sprintf(strerror,"Outmost layer is input as %d, but the surface is not cation terminated.\n", outmost_material);
              nerror (strerror);
            }
            else if ((6 == outmost_material) && (31 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 33; // Outmost: GaAs; Center: Ga. Replace with As. 
            else if ((6 == outmost_material) && (33 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 31; // Outmost: GaAs; Center: As. Replace with Ga. 
            else if ((6 == outmost_material) && (ctr_atom_natyp != 31) && (ctr_atom_natyp != 33)) {
              sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmost_material, ctr_atom_natyp);
              nerror (strerror);
            }
            else if ((7 == outmost_material) && (30 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 34; // ZnSe, Center-Zn, Replace with Se. 
            else if ((7 == outmost_material) && (34 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 30; // ZnSe, Center-Se, Replace with Zn. 
            else if ((7 == outmost_material) && (ctr_atom_natyp != 34) && (ctr_atom_natyp != 30)) {
              sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmost_material, ctr_atom_natyp);
              nerror (strerror);
            }
            else if ((8 == outmost_material) && (30 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 16; // ZnS, Center-Zn, Replace with S. 
            else if ((8 == outmost_material) && (16 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 30; // ZnS, Center-S, Replace with Zn. 
            else if ((8 == outmost_material) && (ctr_atom_natyp != 30) && (ctr_atom_natyp != 16)) {
              sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmost_material, ctr_atom_natyp);
              nerror (strerror);
            }
            else if ((9 == outmost_material) && ((33 == ctr_atom_natyp)||(15 == ctr_atom_natyp))) neighbors_natyp[jneighbor] = 31; // Outmost: alloyGaAsP; Center: As or P. Replace with Ga. 
            else if ((9 == outmost_material) && (31 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 33; // Outmost: alloyGaAsP; Center: Ga. Replace with As. This is a completely random choice. 
            else if ((9 == outmost_material) && (ctr_atom_natyp != 31) && (ctr_atom_natyp != 33) && (ctr_atom_natyp != 15)) {
              sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmost_material, ctr_atom_natyp);
              nerror (strerror);
            }
            else if ((11 == outmost_material) && (31 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 15; // Outmost: GaP; Center: Ga. Replace with P. 
            else if ((11 == outmost_material) && (15 == ctr_atom_natyp)) neighbors_natyp[jneighbor] = 31; // Outmost: GaP; Center: P. Replace with Ga. 
            else if ((11 == outmost_material) && (ctr_atom_natyp != 31) && (ctr_atom_natyp != 15)) {
              sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmost_material, ctr_atom_natyp);
              nerror (strerror);
            }
		      }
          // Chiral ligands should not contribute to strain
          if ((4 == neighbors_natyp[jneighbor]) || (6 == neighbors_natyp[jneighbor]) || (7 == neighbors_natyp[jneighbor])){
            continue;
          } else{
            // Default path!
          bond_lengths[jneighbor] = ret_ideal_bond_len(ctr_atom_natyp, neighbors_natyp[jneighbor], crystal_structure); 
          }
          // mpi_print("ctr_atom_natyp = %d, jneighbor = %d, neighbors_natyp[jneighbor] = %d, bond_lengths[jneighbor]=%g \n", ctr_atom_natyp, jneighbor, neighbors_natyp[jneighbor], bond_lengths[jneighbor]); 
	      }
        
        tetrahedron_vol_ref[jatom] = calc_regular_tetrahedron_volume(bond_lengths[0], bond_lengths[1], bond_lengths[2], bond_lengths[3]);
        // mpi_print("The reference tetrahedron volume is %.7f \n", tetrahedron_vol_ref[jatom]);  
			}
		}
		fclose(pf);
	}
	else {
		fprintf(stderr, "\n\nNo file allNeighborBonds.par in current working directory -- the program is exiting!!!\n\n");
		fflush(stdout);
		exit(EXIT_FAILURE);
	}
  
	return;
}

/****************************************************************************/
/* compute strain scaling factor for each atom
 * strain_scale = 1. + a4 * (Omega/Omega0 - 1),
 * where Omega is volume of tetrahedron formed by atom's nearest neighbors */

void calc_strain_scale(double *strain_scale, vector *atom_neighbors, double *tetrahedron_vol_ref, atom_info *atom, double *a4_params, double *a5_params, long natoms){

	FILE *pf;
	long jatom;
	vector v1, v2, v3;//, v4;
	char ctr_atom[3];
	double tmp2, tmp3, *tetrahedron_vol, strain;
	int tmp1;

	if ((tetrahedron_vol = (double *) calloc(natoms, sizeof(double))) == NULL) nerror("tetrahedron_vol");

	if ( access("strain.par", F_OK) != -1 ) {
	
		//printf("Reading strain scale from strain.par...\n");    
		pf = fopen("strain.par", "r");
		for (jatom = 0; jatom < natoms; jatom++) {
			fscanf(pf, "%i %s %lg %lg %lg %lg %lg", &tmp1, ctr_atom, &tmp2, &tmp3, 
					&tetrahedron_vol[jatom], &tetrahedron_vol_ref[jatom], &strain_scale[jatom]);
		}
		fclose(pf);

	}
	else {
		// TODO: parallelize
		for (jatom = 0; jatom < natoms; jatom++) {

			if ((atom[jatom].Zval == 2) || (atom[jatom].Zval == 3)) { 
				strain_scale[jatom] = 1.0; // Make passivation ligands strainless
			}
			else {
				v1 = retSubtractedVectors(atom_neighbors[4*jatom], atom_neighbors[4*jatom+3]);
				v2 = retSubtractedVectors(atom_neighbors[4*jatom+1], atom_neighbors[4*jatom+3]);
				v3 = retSubtractedVectors(atom_neighbors[4*jatom+2], atom_neighbors[4*jatom+3]);

				tetrahedron_vol[jatom] = fabs(retDotProduct(v1, retCrossProduct(v2, v3)))/6.;
				strain = (tetrahedron_vol[jatom]/tetrahedron_vol_ref[jatom] - 1.);
				strain_scale[jatom] = (1. + a4_params[atom[jatom].idx]*strain + a5_params[atom[jatom].idx]*cube(strain));
			}
		}
	}

	pf = fopen("strain.dat","w");
	for (jatom = 0; jatom < natoms; jatom++) {
		fprintf(pf, "%ld %s %.8f %.8f %.8f %.8f %.8f\n", jatom, atom[jatom].atyp, a4_params[atom[jatom].idx], a5_params[atom[jatom].idx], 
				tetrahedron_vol[jatom], tetrahedron_vol_ref[jatom], strain_scale[jatom]);
	}
	fclose(pf);

	free(tetrahedron_vol);

	return;
}

/****************************************************************************/
double calc_regular_tetrahedron_volume(double bond_length1, double bond_length2, double bond_length3, double bond_length4) {
	
	vector v1, v2, v3, v4; 
	v1.x=sqrt(8./9.);   v1.y=0.;   v1.z=-1./3.;
	v2.x=-sqrt(2./9.);   v2.y=sqrt(2./3.);   v2.z=-1./3.;
	v3.x=-sqrt(2./9.);   v3.y=-sqrt(2./3.);   v3.z=-1./3.;
	v4.x=0.;   v4.y=0.;   v4.z=1.;
	
	v1 = retScaledVector(v1, bond_length1);
	v2 = retScaledVector(v2, bond_length2);
	v3 = retScaledVector(v3, bond_length3);
	v4 = retScaledVector(v4, bond_length4);

	vector side1 = retSubtractedVectors(v1, v4); 
	vector side2 = retSubtractedVectors(v2, v4); 
	vector side3 = retSubtractedVectors(v3, v4); 

	double regular_tetrahedron_volume = fabs(retDotProduct(side1, retCrossProduct(side2, side3)))/6.;

	return regular_tetrahedron_volume; 
}

/****************************************************************************/
