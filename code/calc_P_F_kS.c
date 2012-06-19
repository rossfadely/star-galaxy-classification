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
// calc_P_F_kS.c
//
// Calculate the likelihood of each star for each template, marginalized 
// over the coefficient of the fit. 
//

#include "HBSGsep.h"

void calc_P_F_kS(long ii) {

    //External variables
    extern double P_floor,fluxunitfactor,probfrac;
    extern long Ncstep;
	
    // Internal variables
    long jj,kk,goodflag[Nstartemplate],count=0;
    double h,w,lh,rh;
    double cval,lncval,P_F_CkSmax,P_F_CkS,P_lnC_kS;
	double P_C_kS,temp_P_F_KS[Nstartemplate];
    double coeff,coefferr,chi2,max_P_F_kS=1e-300,min_chi2=1e300;
    double dflux[Nfilter],dfluxerr[Nfilter],mflux[Nfilter],detdfluxerr;
	
	
    // Loop over the templates to calculate P_F_kS for each template
    for (jj=0; jj<Nstartemplate; jj++) {
		
        temp_P_F_KS[jj]=0.0;
		
		// Fit the template using linear techniques
		for (kk=0; kk<Nfilter; kk++) {
				
			mflux[kk]    = modelflux_stars[kk][jj];
			dflux[kk]    = dataflux[kk][ii];
			dfluxerr[kk] = datafluxerr[kk][ii];
				
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
			
		// Calculate chi2 for above fit
		for (kk=0; kk<Nfilter; kk++) {
			
			if (kk==0) 
				chi2 = 0.0;
	    
			chi2 += (dflux[kk] - coeff * mflux[kk]) * \
					(dflux[kk] - coeff * mflux[kk]) / \
					dfluxerr[kk] / dfluxerr[kk];
				
		}
						
		// Integrate over 3sigma of uncertainty of coefficient
		// using h size steps, return the marg. like. of template
		for (kk=0; kk<Ncstep; kk++) {
				
			h = coefferr * 3 / (Ncstep-1) * 2.0 ;
				
			if (kk == 0 || kk == Ncstep-2)
				w=3.0/8.0;
			else if (kk == 1 || kk == Ncstep-3)
				w=7.0/6.0;
			else if (kk == 2 || kk == Ncstep-4)
				w=23.0/24.0;
			else w=1.0;
				
			cval   = coeff - h * (Ncstep-1) / 2.0 + h * kk;
			lncval = log(cval);
				
			if (cval>0 && chi2<1400) {
				P_F_CkSmax = 1.0 / sqrt(2.0 * 3.14159) / detdfluxerr * \
							 exp(-0.5 * chi2);
				P_F_CkS    = P_F_CkSmax * exp(-0.5 * (cval - coeff) * \
											  (cval - coeff) / \
											  coefferr / coefferr);
				
				P_lnC_kS   = 1.0/sqrt(2*3.14159* star_coeff_var[jj]) * \
						     exp(-0.5*(lncval - star_coeff_mean[jj]) * \
								 (lncval - star_coeff_mean[jj]) / \
								 star_coeff_var[jj]);
				P_C_kS     = P_lnC_kS / cval;
		
				temp_P_F_KS[jj] += w * h * P_F_CkS * P_C_kS;	
			} else {
				P_F_CkS = P_floor;
				P_C_kS  = P_floor;
				
				temp_P_F_KS[jj] += P_floor;
			}
		}

		if (temp_P_F_KS[jj]>max_P_F_kS) {
			max_P_F_kS=temp_P_F_KS[jj];
		}
		if (chi2<min_chi2) {
			min_chi2 = chi2;
		}		
		
    }
	
	
    //Alloc memory for number of templates
    if(!(starsparse[0][ii] = (double *)calloc(1,sizeof(double)))) {
		printf("Failed to allocate template sparse memory in calc_P_F_kS\n");
		return;
	}
	
	
    // Loop over the templates to N_P_F_kS and goodflag for each template
    for (jj=0; jj<Nstartemplate; jj++) {
				
        goodflag[jj] = 0;
		
		if (temp_P_F_KS[jj]/max_P_F_kS > probfrac) {
			starsparse[0][ii][0] += 1;
			goodflag[jj] = 1;
		}
		
    }
	
    //Alloc space for the indicies and the P_F_kS values
    if(!(starsparse[1][ii]=(double *)calloc(starsparse[0][ii][0],sizeof(double)))) {
		printf("Failed to allocate index sparse memory in calc_P_F_kS\n");
		return;
	}
    if(!(starsparse[2][ii]=(double *)calloc(starsparse[0][ii][0],sizeof(double)))) {
		printf("Failed to allocate value sparse memory in calc_P_F_kS\n");
		return;
	}
	

    for (jj=0; jj<Nstartemplate; jj++) {
		
        if (goodflag[jj]==1) {
			starsparse[1][ii][count] = (double)jj;
			starsparse[2][ii][count] = temp_P_F_KS[jj];
			count+=1;
		}

    }
	
    star_minchi[ii] = min_chi2;

    if ((ii % 500)==0) {
        printf("\rCalculating star sparse array for the %ldth data point",ii);		
        fflush( stdout );
		sleep(1);
    }
	
	
}
