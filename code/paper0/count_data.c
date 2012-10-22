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
// count_data.c
//
// Count number of sources in input data fits file. 
//

#include "HBSGsep.h"

void count_data(char *datafile,long *N) {
	
    fitsfile *fptr;       
    int status, hdunum = 2, hdutype;
    long longnull, nrow;
    float floatnull;
    char strnull[10]; 
	
    strcpy(strnull, " ");
    longnull  = 0;
    floatnull = 0.;
    status    = 0;
	
    if (fits_open_file(&fptr, datafile, READONLY, &status))
        printfitserror(status);
	
    if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status)) 
        printfitserror(status);
	
    fits_get_num_rows(fptr,&nrow,&status);
	
    *N = nrow;
	
    if (fits_close_file(fptr, &status)) 
        printfitserror(status);
}


