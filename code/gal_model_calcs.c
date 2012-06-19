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
// gal_model_calcs.c 
//
// Read in galaxy SEDs, regrid them, calculate the model flux in each 
// filter for each redshift.
//


#include "HBSGsep.h"

void gal_model_calcs(void) {
	
    //External variables
    extern char galssedinput[FILEPATH_LENGTH];
    extern double zmin,zmax;
    extern long Nz;
	
    //Internal variables
    long ii,jj,kk;
    long sedlength;
    long *N;
    double z;
    double *plam,*pval,*pvalfine,*modval;

    N = (long *)malloc(sizeof(long));
	
    //How many gal SEDS?
    printf("\nUsing %ld gal templates over %ld redshifts\n", \
		   Ngaltemplate,Nz);
	
    modval=(double *)malloc(sizeof(double));
    modelflux_gals=malloc(Nfilter*sizeof(double*));
	
	
    // Calculate Model fluxes
    zstep=(zmax-zmin)/(Nz-1);
		
    for (ii=0;ii<Ngaltemplate;ii++) {	

        get_filelength(ii,galssedinput, N);
        sedlength = *N;

		plam=(double *)malloc(*N * sizeof(double));
		pval=(double *)malloc(*N * sizeof(double));
	
		read_file(ii,galssedinput,plam,pval);
	
		for (jj=0; jj<Nz; jj++) {
			z=zmin+jj*zstep;
			for (kk=0; kk<Nfilter; kk++) {
				if (ii==0 && jj==0) {
					modelflux_gals[kk] = (double *)malloc(Ngaltemplate * Nz * \
														  sizeof(double));
				}
				pvalfine = (double *)malloc(filter_lgth_fine[kk] * sizeof(double));
				
				regrid_sed(z,plam,pval,filter_lgth_fine[kk],sedlength, \
						   filter_lamb_fine[kk],pvalfine);
			
				integrate_sed(filter_lgth_fine[kk],pvalfine,filter_lamb_fine[kk],\
							  filter_thru_fine[kk],modval);
			
				modelflux_gals[kk][jj+ii*Nz]=*modval;
				
				free(pvalfine);
			}
		}
		
		free(plam);
		free(pval);
    }
    free(modval);
	
}


