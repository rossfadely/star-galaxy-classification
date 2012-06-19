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
//  coeff_calcs.c
//
//  Fit star and galaxy templates, determine the coefficient priors.
//

#include "HBSGsep.h"

void coeff_calcs(void) {
	
    //External variables
    extern long Nz;
	
    //Internal variables
    long ii,jj;

    printf("\nCalculating values for coefficient priors...\n\n ");
	
    double *star_chi2;
    double *star_lncoeff,*star_coeffvar;
 	
    star_coeff_mean = (double *)malloc(Nstartemplate*sizeof(double));
    star_coeff_var  = (double *)malloc(Nstartemplate*sizeof(double));
    star_chi2       = (double *)malloc(Ndata*sizeof(double));
    star_lncoeff    = (double *)malloc(Ndata*sizeof(double));
    star_coeffvar   = (double *)malloc(Ndata*sizeof(double));
	
    for (ii=0;ii<Nstartemplate;ii++) {
		
        fit_star_template(ii,star_lncoeff,star_coeffvar,star_chi2);
        coeff_meanvar(1,ii,star_lncoeff,star_coeffvar,star_chi2);
		
    }
    free(star_lncoeff);
    free(star_coeffvar);
    free(star_chi2);
		
    double *gal_chi2;
    double *gal_lncoeff,*gal_coeffvar;
	
    gal_coeff_mean = (double *)malloc(Ngaltemplate*Nz*sizeof(double));
    gal_coeff_var  = (double *)malloc(Ngaltemplate*Nz*sizeof(double));
    gal_chi2       = (double *)malloc(Ndata*sizeof(double));
    gal_lncoeff    = (double *)malloc(Ndata*sizeof(double));
    gal_coeffvar   = (double *)malloc(Ndata*sizeof(double));
	
	
    for (ii=0;ii<Ngaltemplate;ii++) {
        for (jj=0; jj<Nz; jj++) {
			
            fit_gal_template(ii,jj,gal_lncoeff,gal_coeffvar,gal_chi2);
            coeff_meanvar(0,ii*Nz+jj,gal_lncoeff,gal_coeffvar,gal_chi2);
			
        }
    }
    free(gal_lncoeff);
    free(gal_coeffvar);
    free(gal_chi2);

}


