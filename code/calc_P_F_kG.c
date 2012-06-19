//
// This is part of HBSGsep, a hierarchical routine to classify stars and 
// galaxies using photometric data.
//
// Please consult http://github.com/rossfadely/star-galaxy-classifiction
// for associated documentation.
//
//
// ======================================================================
//
// calc_P_F_kG.c
//
// Calculate the likelihood of each star for each template, marginalized 
// over the coefficient of the fit and redshift. This assumes a flat 
// redshift prior.
//


#include "HBSGsep.h"

void calc_P_F_kG(long ii) {
	
	
    //External variables
    extern double P_floor,fluxunitfactor,probfrac;
    extern double zmin,zstep;
    extern long Nz,Ncstep;
	
    // Internal variables
    long jj,ll,kk,count=0;
    double h,w,lh,rh,z;
    double coeff,coefferr,chi2,min_chi2=1e300;
    double dflux[Nfilter],dfluxerr[Nfilter],mflux[Nfilter];
    double detdfluxerr,max_P_F_kG=1e-300;
    double P_F_CkzGmax=0.0,P_F_CkzG=0.0,P_F_kzG=0.0,P_lnC_kzG=0.0;
    double P_C_kzG=0.0,norm_P_z_G=0.0,temp_P_F_kG[Ngaltemplate];
    double cval,lncval;
    int goodflag[Ngaltemplate];
	
    for (jj=0;jj<Ngaltemplate;jj++) {

        temp_P_F_kG[jj]=0.0;
			
        norm_P_z_G=0.0;
        for (ll=0; ll<Nz; ll++) {

            P_F_kzG = 0.0;

            z = zmin+ll*zstep;
			
            for (kk=0; kk<Nfilter; kk++) {
					
                mflux[kk]    = modelflux_gals[kk][ll+jj*Nz];
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
            // using h size steps, return the marg. like. of template at z
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
                    P_F_CkzGmax = 1.0 / sqrt(2.0 * 3.14159) / detdfluxerr * \
                                  exp(-0.5 * chi2);
                    P_F_CkzG    = P_F_CkzGmax * exp(-0.5 * (cval - coeff) * \
                                                    (cval - coeff) / \
                                                    coefferr / coefferr);
                    P_lnC_kzG   = 1.0/sqrt(2 * 3.14159 * gal_coeff_var[jj*Nz+ll]) * \
                                  exp(-0.5*(lncval - gal_coeff_mean[jj*Nz+ll]) * \
                                      (lncval - gal_coeff_mean[jj*Nz+ll]) / \
                                      gal_coeff_var[jj*Nz+ll]);
                    P_C_kzG     = P_lnC_kzG / cval;
                    P_F_kzG    += w * h * P_F_CkzG * P_C_kzG;	
                } else {
                    P_F_CkzG    = P_floor;
                    P_C_kzG     = P_floor;
                    P_F_kzG    += P_floor;
                }
					
            }
								
            norm_P_z_G      += 1.0 / Nz;
            temp_P_F_kG[jj] += P_F_kzG / Nz;
            
				
            if (chi2<min_chi2) {
                min_chi2=chi2;
            }
		
        }
			
        temp_P_F_kG[jj] /= norm_P_z_G;
			
			
        if (temp_P_F_kG[jj]>max_P_F_kG) {
            max_P_F_kG = temp_P_F_kG[jj];
        }
			
			
    }
		
		
	
    // Loop over the templates and zpriors to N_P_F_kG and goodflag for each template
    if(!(galsparse[0][ii]=(double *)calloc(1,sizeof(double)))) {
        printf("Failed to allocate template sparse memory in calc_P_F_kG\n");
        return;
    }	
	
		
    for (jj=0;jj<Ngaltemplate;jj++) {
        goodflag[jj] = 0;
        if (temp_P_F_kG[jj]/max_P_F_kG > probfrac) {
            galsparse[0][ii][0] += 1;
            goodflag[jj] = 1;
        }	
    }
		

    //Alloc space for the indicies and the P_F_kG values
    if(!(galsparse[1][ii]=(double *)calloc(galsparse[0][ii][0],sizeof(double)))) {
        printf("Failed to allocate index sparse memory in calc_P_F_kG\n");
        return;
    }	
    if(!(galsparse[2][ii]=(double *)calloc(galsparse[0][ii][0],sizeof(double)))) {
        printf("Failed to allocate value sparse memory in calc_P_F_kG\n");
        return;
    }

    count=0;
		
    for (jj=0;jj<Ngaltemplate;jj++) {
        if (goodflag[jj]==1) {
            galsparse[1][ii][count] = jj;
            galsparse[2][ii][count] = temp_P_F_kG[jj];
            count += 1;
        }			
    }
		
	
    gal_minchi[ii] = min_chi2;
	
    if ((ii % 500)==0) {
        printf("\rCalculating galaxy sparse array for the %ldth data point",ii);		
        fflush( stdout );
        sleep(1);
    }
	
}
