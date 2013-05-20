//
// This is part of HBSGsep, a hierarchical routine to classify astronomical 
// sources using photometric data.
//
// Please consult http://github.com/rossfadely/star-galaxy-classifiction
// for associated documentation.
//
//
// =====================================================================
//
//  fit_one.c
//
//  Fit a given source with a single template, record associated coefficient, 
//  the coefficient error, and chi2 for the fit.
//

#include <stdio.h>

int fit_one(long,double*,double*,double*,double*,double*,double*);

void fit_models(long Ndata,long Nmodel,long Nfilter,
                double **models,double **data,double **dataerr,
                double **coeffs,double **coefferrs,double **chi2s) {
    
    int err_flag;
    long ii,jj;
    
    for (ii=0; ii<Ndata; ii++) {
        for (jj=0; jj<Nmodel; jj++) {
            err_flag = fit_one(Nfilter,models[jj],data[ii],dataerr[ii],&coeffs[ii][jj],
                               &coefferrs[ii][jj],&chi2s[ii][jj]);
        }
    }
    
}


int fit_one(long Nfilter,double *modelflux,double *dataflux,
            double *datafluxerr,double *coeff,double *coefferr,
            double *chi2) {
	
    int flag=0;
    long ii;
    double lh=0.0,rh=0.0;
    
    *chi2 = 0.0;
		
    // Fit the model using linear techniques
    for (ii=0; ii<Nfilter; ii++) {
                
        lh +=  modelflux[ii] * modelflux[ii] / datafluxerr[ii] / datafluxerr[ii];
        rh +=  dataflux[ii] * modelflux[ii] / datafluxerr[ii] / datafluxerr[ii];

    }
		
    *coeff = rh/lh;
    *coefferr = sqrt(1.0/lh);
        
    // Calculate chi2 for above fit
    if (*coeff<=0) {
        flag = 1;
    } else {
        for (ii=0; ii<Nfilter; ii++) {
			            
            *chi2 += (dataflux[ii] - *coeff * modelflux[ii]) * \
                     (dataflux[ii] - *coeff * modelflux[ii]) / \
                     datafluxerr[ii] / datafluxerr[ii];
        }
        
    }

    return flag;
}
