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
//  star_model_mags.c
//
//  Calculate model magnitudes of star SEDs.
//

#include "HBSGsep.h"

void star_model_mags(void) {
	
    //External variables
    extern char starmodmagsfile[FILEPATH_LENGTH];
	
    //Internal variables
    long jj,kk;
    double **modmags;
	
    modmags=(double **)malloc(Nfilter * sizeof(double *));
	
    for (kk=0;kk<Nfilter;kk++) {
        modmags[kk] = (double *)malloc(Nstartemplate * sizeof(double));
		for (jj=0; jj<Nstartemplate; jj++) {
			modmags[kk][jj] = (-2.5) * log10(modelflux_stars[kk][jj] / \
											 norm[kk]);
		}
    }
	
    write_modelmags(starmodmagsfile,Nstartemplate,modmags);
		
    for (kk=0;kk<Nfilter;kk++) {
        free(modmags[kk]);
    }
    free(modmags);
}


