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
//  read_hyppars.c
//
//  Read in hyperparameters to initialize optimization.
//
//

#include "HBSGsep.h"

void read_hyppars(char *hypinfile,double *hypparms) {
	
    fitsfile *fptr;       
    int status, hdunum=2, hdutype, anynull;
    long nrow_l;
    float floatnull;
	
    floatnull = 0.;
    status    = 0;
	
    if (fits_open_file(&fptr, hypinfile, READONLY, &status))
        printfitserror(status);
	
    if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status)) 
        printfitserror(status);
	
    fits_get_num_rows(fptr,&nrow_l,&status);
	
    fits_read_col(fptr, TDOUBLE, 1, 1, 1, nrow_l, &floatnull, \
                  hypparms, &anynull, &status);
	
    if (fits_close_file(fptr, &status)) 
        printfitserror(status);
}
