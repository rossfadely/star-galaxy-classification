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
//  fit_star_template.c
//
//  Fit a given star template to all the data, and record the 
//  coefficients, their uncertainties, and the associated chi2.
//

#include "HBSGsep.h"

void fit_star_template(long ii,double *star_lncoeff, \
                       double *star_coeffvar,double *star_chi2) {
	
    //External variables
    extern double fluxunitfactor;
	
    // Internal variables
    long jj,kk;
    double lh,rh;
    double coeff,coefferr,chi2;
    double dflux[Nfilter],dfluxerr[Nfilter],mflux[Nfilter],detdfluxerr;

    // Loop over the data
    for (jj=0; jj<Ndata; jj++) {
		
        // Fit the template using linear techniques
        for (kk=0; kk<Nfilter; kk++) {
				
            mflux[kk]    = modelflux_stars[kk][ii];
            dflux[kk]    = dataflux[kk][jj];
            dfluxerr[kk] = datafluxerr[kk][jj];
				
            dflux[kk]    = dflux[kk] * fluxunitfactor;
            dfluxerr[kk] = dfluxerr[kk] * fluxunitfactor;
				
            if (kk == 0) {
                lh = 0.0;
                rh = 0.0;
                detdfluxerr = dfluxerr[kk];
            } else {
                detdfluxerr *= dfluxerr[kk];
            }
				
            lh +=  mflux[kk] * mflux[kk] / dfluxerr[kk] / dfluxerr[kk];
            rh +=  dflux[kk] * mflux[kk] / dfluxerr[kk] / dfluxerr[kk];
				
            if (kk == Nfilter-1) {
                coeff = rh/lh;
                coefferr = sqrt(1.0/lh);
            }
			
        }
		
        if (coeff<=0)
            coeff = 1.0;
		
        // Calculate chi2 for above fit
        for (kk=0; kk<Nfilter; kk++) {
			
            if (kk==0) 
                chi2 = 0.0;
	    
            chi2 += (dflux[kk] - coeff * mflux[kk]) * \
                    (dflux[kk] - coeff * mflux[kk]) / \
                    dfluxerr[kk] / dfluxerr[kk];
			
        }
		
    star_lncoeff[jj]  = log(coeff);
    star_coeffvar[jj] = coefferr * coefferr;
    star_chi2[jj]     = chi2;
       
    }
}
