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
// star_model_calcs.c 
//
// Read in star SEDs, regrid them, calculate the model flux in each 
// filter for each redshift.
//

#include "HBSGsep.h"

void star_model_calcs(void) {
	
    //External variables
    extern char starssedinput[FILEPATH_LENGTH];
	
    //Internal variables
    long ii,kk;
    long sedlength;
    long *N;
    double *plam,*pval,*pvalfine,*modval;

    N = (long *)malloc(sizeof(long));
	
    //How many star SEDS?
    printf("\nUsing %ld star templates\n",Nstartemplate);
	
    modval          = (double *)malloc(sizeof(double));
    modelflux_stars = malloc(Nfilter*sizeof(double*));
	
    // Calculate Model fluxes
    for (ii=0;ii<Nstartemplate;ii++) {	

        get_filelength(ii,starssedinput,N);
        sedlength = *N;

        plam = (double *)malloc(sedlength * sizeof(double));
        pval = (double *)malloc(sedlength * sizeof(double));

        read_file(ii,starssedinput,plam,pval);

        for (kk=0; kk<Nfilter; kk++) {
            if (ii==0) {
                modelflux_stars[kk] = (double *)malloc(Nstartemplate * \
                                                       sizeof(double));
            }
            pvalfine = (double *)malloc(filter_lgth_fine[kk] * sizeof(double));
			
            regrid_sed(0.0,plam,pval,filter_lgth_fine[kk],sedlength, \
                       filter_lamb_fine[kk],pvalfine);
			
            integrate_sed(filter_lgth_fine[kk],pvalfine,filter_lamb_fine[kk], \
                          filter_thru_fine[kk],modval);
			
            modelflux_stars[kk][ii] = *modval;
			
            free(pvalfine);
        }
        free(plam);
        free(pval);
    }
    free(modval);
	
}


