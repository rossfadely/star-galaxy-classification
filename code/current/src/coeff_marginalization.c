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
//  coeff_marginalization.c
//
//  Marginalize over uncertainty in the fit coefficient.
//

#include <math.h>

double coeff_marginalization(long Nstep,double Nsigma,
                             double coeff,double coefferr,
                             double prior_mean,double prior_var,
                             double detdfluxerr,double chi2){
    
    long ii;
    double h,w;
    double cval,lncval;
    double dNstep,A,B;
    double marg_like=0.0,max_like=0.0,like=0.0;
    double prior=0.0,log_norm_coeff=0.0;
    
    dNstep = (double)Nstep;
	
    A = 1.0 / sqrt(2.0 * 3.14159265);
    B = A / sqrt(prior_var);
    h = 2. *  Nsigma * coefferr  / (dNstep-1.0);
    cval0 = coeff - Nsigma * coefferr;
    max_like = A / detdfluxerr * exp(-0.5 * chi2);
    
    for (ii=0; ii<Nstep; ii++) {
    
        if (ii == 0 || ii == Nstep-2)
            w=3.0/8.0;
        else if (ii == 1 || ii == Nstep-3)
            w=7.0/6.0;
        else if (ii == 2 || ii == Nstep-4)
            w=23.0/24.0;
        else w=1.0;
        
        cval   = cval0 + h * (double)ii;
        lncval = log(cval);
        
        like = max_like * exp(-0.5 * (cval - coeff) * \
                              (cval - coeff) / coefferr / coefferr);
        log_norm_coeff = B * exp(-0.5*(lncval - prior_mean) * \
                         (lncval - prior_mean) / prior_var);
        prior = log_norm_coeff / cval;
        marg_like += h * w * like * prior;
    
    }

    return marg_like;
    
}
