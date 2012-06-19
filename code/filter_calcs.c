//
// This is part of HBSGsep, a hierarchical routine to classify stars and 
// galaxies using photometric data.
//
// Please consult http://github.com/rossfadely/star-galaxy-classifiction
// for associated documentation.
//
//
// ==========================================================================
//
//  filter.c
//
//  Routine reads in the model seds and filter throughput response curves,
//  regrids the filters onto a finer grid for interpolation, and calls a 
//  routine to calculate the filter flux zeropoints.  Testing indicates 
//  details of regridding and interpolation scheme is not very important.
//

#include "HBSGsep.h"

void filter_calcs(void) {
	
	//External variables
	extern char filtersinput[FILEPATH_LENGTH];
	extern char galssedinput[FILEPATH_LENGTH];
	extern char starssedinput[FILEPATH_LENGTH];
	
	//Internal variables
	long *N,ii;
	long sedlength = 0;
	long regridfactor;
	long filterlength[Nfilter]; 

	double *pnorm;
	double *filtlamb,*filtthru;

	//Allocate temporary variables
	N     = (long *)malloc(sizeof(long));
	pnorm = (double *)malloc(sizeof(double));
	
	//How many filters?
	get_num_files(filtersinput, N);
	Nfilter = *N;
	if (Nfilter<2) {
	    printf("Need more than 1 filter\n");
		return;
	}
	printf("\n\nFound %ld Filters\n",Nfilter);
	
	//Read in the number of star/galaxy SEDs
	get_num_files(starssedinput, N);
	Nstartemplate = *N;
	get_num_files(galssedinput, N);
	Ngaltemplate  = *N;
	
	//Find the finest SED amongst the bunch
	for (ii=0;ii<Nstartemplate;ii++) {	
		get_filelength(ii,starssedinput,N);
		if (*N * 2 > sedlength) 
			sedlength = *N * 2;
	}
	for (ii=0;ii<Ngaltemplate;ii++) {	
		get_filelength(ii,galssedinput,N);
		if (*N * 2 > sedlength)
			sedlength = *N * 2;
	}
	
	
	//Allocate final filter arrays, which are globals
	filter_lgth_fine = (long *)malloc(Nfilter*sizeof(long)); 
	filter_lamb_fine = malloc(Nfilter*sizeof(double*));
	filter_thru_fine = malloc(Nfilter*sizeof(double*));
	norm             = (double *)malloc(Nfilter*sizeof(double));
	
	
	//Loop over the filters in the file
	for (ii=0;ii<Nfilter;ii++) {
		
		//get length
		get_filelength(ii,filtersinput, N);
		filterlength[ii]     = *N;
		regridfactor = round((float)sedlength / (float)*N);
		filter_lgth_fine[ii] = *N * regridfactor;
		
		//alloc filter arrays
		filtlamb = (double *)malloc(*N * sizeof(double));
		filtthru = (double *)malloc(*N * sizeof(double));
		filter_lamb_fine[ii] = (double *)malloc(regridfactor * *N * sizeof(double));
		filter_thru_fine[ii] = (double *)malloc(regridfactor * *N * sizeof(double));
		
		//read in the 2 column ascii filter file
		read_file(ii,filtersinput,filtlamb,filtthru);
		
		//regrid the filter to user spec, using gsl spline interpolation
		regrid_filter(filtlamb,filtthru,filterlength[ii],filter_lgth_fine[ii], \
					  filter_lamb_fine[ii],filter_thru_fine[ii]);
		
		//calculate the flux zeropoint
		calc_normalization(filter_lamb_fine[ii],filter_thru_fine[ii], \
						   filter_lgth_fine[ii],pnorm);
		norm[ii] = *pnorm;
		printf("Filter %ld has (AB) zeropoint flux normalization: %g\n",ii,norm[ii]);
		
		free(filtlamb);
		free(filtthru);
	}
	free(pnorm);
	free(N);
	

}


