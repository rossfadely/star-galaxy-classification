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
//  gal_model_mags.c
//
//  Calculate model magnitudes of galaxy SEDs.
//

#include "HBSGsep.h"

void gal_model_mags(void) {
	
    //External variables
    extern long Nz;
    extern char galmodmagsfile[FILEPATH_LENGTH];
	
    //Internal variables
    long jj,kk,ll;
    double **modmags;
	
    modmags=(double **)malloc(Nfilter * sizeof(double *));
	
    for (kk=0;kk<Nfilter;kk++) {
        modmags[kk]=(double *)malloc(Ngaltemplate * Nz * sizeof(double));
		for (jj=0; jj<Ngaltemplate; jj++) {
			for (ll=0; ll<Nz; ll++) {
				modmags[kk][ll+jj*Nz] = (-2.5) * \
										log10(modelflux_gals[kk][ll+jj*Nz] / \
											  norm[kk]);
			}
		}
    }
	
    write_modelmags(galmodmagsfile,Ngaltemplate * Nz,modmags);
		
    for (kk=0;kk<Nfilter;kk++) {
      free(modmags[kk]);
    }
    free(modmags);
}


