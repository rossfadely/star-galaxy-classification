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
// calc_datavals.c
//
//  Read in the FITS data.  FITS file should have magnitudes in the 
//  first Nfilter columns, then the uncertainties, then extra columns 
//  (not read in).  Data and uncertainties are converted to flux units, 
//  with special treatment for missing or undetected data.  The latter 
//  needs to be taylored to the catalog, or vice versa.
//

#include "HBSGsep.h"

void calc_datavals(long ii,char *datafile,double norm,double *flux, \
				   double *fluxerr) {
	
    //external variables
    extern double noisefudge;
	
    //internal variables
    fitsfile *fptr;       
    int status, hdunum=2, hdutype, anynull, jj;
    long longnull, nrow;
    float floatnull;
    char strnull[10]; 
    double *pmag,*perr;
	
	// initialize values for fits files, open and read them.
    strcpy(strnull, " ");
    longnull  = 0;
    floatnull = 0.;
    status    = 0;
	
    if (fits_open_file(&fptr, datafile, READONLY, &status))
        printfitserror(status);
	
	if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status)) 
	    printfitserror(status);
	
	fits_get_num_rows(fptr,&nrow,&status);

	pmag = (double *)malloc(nrow * sizeof(double));
	perr = (double *)malloc(nrow * sizeof(double));
	
	fits_read_col(fptr, TDOUBLE, ii+1, 1, 1, nrow, &floatnull, pmag, \
				  &anynull, &status);
	fits_read_col(fptr, TDOUBLE, ii+Nfilter+1, 1, 1, nrow, &floatnull, \
				  perr, &anynull, &status);

	// Take input magnitudes and errors and convert to flux
	for (jj=0; jj<nrow; jj++) {
	    // Missing data   
	    if (pmag[jj] < 0) {
	        flux[jj]    = pow(10.0,-0.4 * 10.0) * norm;
			fluxerr[jj] = pow(10.0,-0.4 * 10.0) * norm * log(10.0) * \
					      0.4 * 1000000.;
	    // Undetected 
	    } else if (pmag[jj] > 50) {
	        flux[jj]    = 0.0;
			fluxerr[jj] = pow(10.0,-0.4 * perr[jj]) * norm * 2.0;
	    // All ok
	    } else {
	        flux[jj]    = pow(10.0,-0.4 * pmag[jj]) * norm;
			fluxerr[jj] = pow(10.0,-0.4 * pmag[jj]) * norm * log(10.0) \
			              * 0.4 * perr[jj];
	    }
	    // Add noise model
	    fluxerr[jj] = sqrt(pow(fluxerr[jj],2) + \
						   pow(flux[jj] * noisefudge,2));
	}

    if (fits_close_file(fptr, &status)) 
        printfitserror(status);
	
    free(pmag);
    free(perr);
}


