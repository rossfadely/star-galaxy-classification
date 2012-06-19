//
// This is part of HBSGsep, a hierarchical routine to classify stars and 
// galaxies using photometric data.
//
// Please consult http://github.com/rossfadely/star-galaxy-classifiction
// for associated documentation.
//
//
// ==========================================================================
//
// calc_P_F_G.c
//
// Calculate the likelihood of each object to be any type of galaxy 
// template, weighted by the hyperparameters.
//
//

#include "HBSGsep.h"

void calc_P_F_G(long ii,double *hypparms_gals,double *P_F_G) {
	
	
	long jj,index;
	
	for (jj=0; jj<(long)galsparse[0][ii][0]; jj++) {
		
		index = (long)galsparse[1][ii][jj];
	
		P_F_G[ii] += hypparms_gals[index] * galsparse[2][ii][jj]; 

	}
	
}
