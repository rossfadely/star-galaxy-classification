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
// coeff_meanvar.c
//
// Calculate the mean and variance of coefficients produced by fitting 
// templates to all the data.  
//

#include "HBSGsep.h"

void coeff_meanvar(int flag,long ii,double *lncoeff_in,double *coeffvar, \
				   double *chi2) {

    long jj;
    double *weight,*lncoeff;
    weight  = (double *)malloc(Ndata*sizeof(double));
    lncoeff = (double *)malloc(Ndata*sizeof(double));

    // Loop over the data
    for (jj=0; jj<Ndata; jj++) {
        weight[jj]  = 1.0 / coeffvar[jj];
		lncoeff[jj] = lncoeff_in[jj];
    }
	
    if (flag==0) {
        gal_coeff_mean[ii] = gsl_stats_wmean(weight,1,lncoeff,1,Ndata);
		gal_coeff_var[ii]  = gsl_stats_wvariance(weight,1,lncoeff,1,Ndata);
		
		if (Ndata<10)
			gal_coeff_var[ii]  = gal_coeff_mean[ii] * -0.5;
		
		if (gal_coeff_var[ii]<0.0)
			gal_coeff_var[ii] *= -1.0;
    } 
    if (flag==1) {
        star_coeff_mean[ii] = gsl_stats_wmean(weight,1,lncoeff,1,Ndata);
		star_coeff_var[ii]  = gsl_stats_wvariance(weight,1,lncoeff,1,Ndata);
		
		if (Ndata<10)
			star_coeff_var[ii]  = star_coeff_mean[ii] * -0.5;
		
		if (star_coeff_var[ii]<0.0)
			star_coeff_var[ii] *= -1.0;
    }
	
    free(weight);
    free(lncoeff);
}
