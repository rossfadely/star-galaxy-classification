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

#include <stdio.h>
#include <math.h>

void single_coeff_marginalization(long,double,double,double,double,
                                  double,double,double,double*);

void coeff_marginalization(long Nstep,long Ndata,
                           long Nmodel,double Nsigma,
                           double *minchi2s,double Dchi2cut,
                           double **coeffs,double **coefferrs,
                           double *prior_means,double *prior_vars,
                           double *detdfluxerr,double **chi2s,
                           double **marg_likes) {

    long ii,jj;
    for (ii=0; ii<Ndata; ii++) {
        for (jj=0; jj<Nmodel; jj++) {
            if ((chi2s[ii][jj]-minchi2s[ii])<Dchi2cut) {
                single_coeff_marginalization(Nstep,Nsigma,coeffs[ii][jj],
                                             coefferrs[ii][jj],prior_means[jj],
                                             prior_vars[jj],detdfluxerr[ii],
                                             chi2s[ii][jj],&marg_likes[ii][jj]);
            }
        }
    }
    
}


void single_coeff_marginalization(long Nstep,double Nsigma,
                                  double coeff,double coefferr,
                                  double prior_mean,double prior_var,
                                  double detdfluxerr,double chi2,
                                  double *marg_like){
    
    long ii;
    double h,w;
    double cval,cval0,lncval;
    double dNstep,A,B;
    double max_like=0.0,like=0.0;
    double prior=0.0,log_norm_coeff=0.0;
    
    dNstep = (double)Nstep;
    *marg_like = 0.0;
	
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
        
        if (cval>0) {
            like = max_like * exp(-0.5 * (cval - coeff) * \
                                  (cval - coeff) / coefferr / coefferr);
            log_norm_coeff = B * exp(-0.5*(lncval - prior_mean) * \
                                     (lncval - prior_mean) / prior_var);
            prior = log_norm_coeff / cval;
            *marg_like += h * w * like * prior;            
        }
    
    }
    
}
