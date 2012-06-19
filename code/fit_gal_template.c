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
//  fit_galaxy_template.c
//
//  Fit a given galaxy template to all the data, and record the 
//  coefficients, their uncertainties, and the associated chi2.
//

#include "HBSGsep.h"

void fit_gal_template(long ii,long jj,double *gal_lncoeff, \
					  double *gal_coeffvar,double *gal_chi2) {
	
	
    //External variables
    extern long Nz;
    extern double fluxunitfactor;
	
    // Internal variables
    long   ll,kk;
    double lh,rh;
    double coeff,coefferr,chi2;
    double dflux[Nfilter],dfluxerr[Nfilter],mflux[Nfilter],detdfluxerr;
		
    for (ll=0; ll<Ndata; ll++) {
        // Calculate best fit
        for (kk=0; kk<Nfilter; kk++) {
			
	    mflux[kk]    = modelflux_gals[kk][jj+ii*Nz];
	    dflux[kk]    = dataflux[kk][ll];
	    dfluxerr[kk] = datafluxerr[kk][ll];
						
	    dflux[kk]    = dflux[kk]*fluxunitfactor;
	    dfluxerr[kk] = dfluxerr[kk]*fluxunitfactor;
			
	    if (kk == 0) {
	        lh=0.0;
		rh=0.0;
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
		
	gal_lncoeff[ll]  = log(coeff);
	gal_coeffvar[ll] = coefferr * coefferr;
	gal_chi2[ll]     = chi2;
				
    }
}
