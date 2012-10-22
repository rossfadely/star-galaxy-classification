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
// printfitserror.c
//
// CFITSIO error messages.
//


#include "HBSGsep.h"

void printfitserror(int status) {

    if (status){
        fits_report_error(stderr, status); 	
        exit(status);
    }
    return;
}
