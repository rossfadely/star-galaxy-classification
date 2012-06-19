//
// This is part of HBSGsep, a hierarchical routine to classify stars and 
// galaxies using photometric data.
//
// Please consult http://github.com/rossfadely/star-galaxy-classifiction
// for associated documentation.
//
//
// =====================================================================
//
// sparce_calcs.c
//
// Call functions to save the marginalized likelihood (over coefficient)
// for all templates above the probability threshold.
// 

#include "HBSGsep.h"

void sparse_calcs(void) {
	
    long ii;
	
    starsparse = (double ***)malloc(3*sizeof(double **));
    for (ii=0; ii<3; ii++) {
        if(!(starsparse[ii] = (double **)malloc(Ndata*sizeof(double *)))) {
            printf("Failed to allocate star sparse memory in sparse_calcs\n");
            return;
        }
    }
	
	
    for (ii=0; ii<Ndata; ii++) {
        calc_P_F_kS(ii);
    }

    printf("\n\n");
	
    galsparse=(double ***)malloc(3*sizeof(double **));
    for (ii=0; ii<3; ii++) {
        if(!(galsparse[ii] = (double **)malloc(Ndata*sizeof(double *)))) {
            printf("Failed to allocate galaxy sparse memory in sparse_calcs\n");
            return;
        }
    }
	
    for (ii=0; ii<Ndata; ii++) {
        calc_P_F_kG(ii);
    }
	
}


