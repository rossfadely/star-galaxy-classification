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
// data_calcs.c
//
// Routine to assign data arrays, and call calc_datavals to calculate 
// flux values.
//

#include "HBSGsep.h"

void data_calcs(void) {
	
	//External variables
	extern char datafile[FILEPATH_LENGTH];
	
	//Internal variables
	long ii,jj,*N;
	double *flux, *fluxerr;
	N = (long *)malloc(sizeof(long));
	
	//Read data mags and uncertainties
	count_data(datafile,N);
	Ndata = *N;
	printf("\nThere are %ld sources in data file\n",Ndata);
	
	//alloc
	dataflux    = malloc(Nfilter*sizeof(double*));
	datafluxerr = malloc(Nfilter*sizeof(double*));
	flux        = (double *)malloc(Ndata * sizeof(double));
	fluxerr     = (double *)malloc(Ndata * sizeof(double));
	
	//read in data for each filter
	for (ii=0;ii<Nfilter;ii++) {

		//alloc 
		dataflux[ii]    = (double *)malloc(Ndata * sizeof(double));
		datafluxerr[ii] = (double *)malloc(Ndata * sizeof(double));
		
		//read in data for ii-th filter
		calc_datavals(ii,datafile,norm[ii],flux,fluxerr);
		for (jj=0; jj<Ndata; jj++) {
			dataflux[ii][jj]    = flux[jj];
			datafluxerr[ii][jj] = fluxerr[jj];
		}		
	}	
	free(N);
	free(flux);
	free(fluxerr);
	
}


