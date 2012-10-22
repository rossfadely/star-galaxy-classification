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
// calc_P_F_S.c
//
// Calculate the likelihood of each object to be any type of star 
// template, weighted by the hyperparameters.
//
//


#include "HBSGsep.h"

void calc_P_F_S(long ii,double *hypparms_stars,double *P_F_S) {
	
	
    long jj,index;
    
    for (jj=0; jj<(long)starsparse[0][ii][0]; jj++) {
		
        index = (long)starsparse[1][ii][jj];
		
        P_F_S[ii] += hypparms_stars[index] * starsparse[2][ii][jj]; 

    }
	
}
