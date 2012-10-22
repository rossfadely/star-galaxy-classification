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
//  fit_calcs.c
//
//  Call routines to fit templates to the data, and record likelihood 
//  related info.
//

#include "HBSGsep.h"

void fit_calcs(void) {
	
    coeff_calcs();
	
    star_minchi = malloc(Ndata*sizeof(double));
    gal_minchi  = malloc(Ndata*sizeof(double));

    sparse_calcs();
	
}


