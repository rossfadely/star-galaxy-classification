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
// write_hyperparms.c
//
// Write hyperparameters to file.
//
//

#include "HBSGsep.h"

void write_hyppars(char *hypoutfile,double *hypparms) {
		
	
    int bitpix = USHORT_IMG; 
    long naxis = 1;
    long naxes[1] = {1};
    int status=0,tfields=1,hdutype;
    char *ttype[1],*tform[1],*tunit[1];          
    char extname[] = " ";           
    fitsfile *fptr;       
	
	
    ttype[0] = (char *)malloc(10 * sizeof(char));
    tform[0] = (char *)malloc(10 * sizeof(char));
    tunit[0] = (char *)malloc(10 * sizeof(char));
	
    strcpy(*(ttype+0),"hyp");
    strcpy(*(tform+0),"1D");
    strcpy(*(tunit+0),"      ");
	
    remove(hypoutfile);   
	
    if (fits_create_file(&fptr,hypoutfile,&status)) 
        printfitserror(status );           
	
    if (fits_create_img(fptr,bitpix,naxis,naxes,&status))
        printfitserror(status);   
	
    if (fits_movabs_hdu(fptr,1,&hdutype,&status))
        printfitserror(status);
	
    if (fits_create_tbl(fptr,BINARY_TBL,Nhyperparms,tfields,ttype,tform,tunit, \
                        extname,&status))
        printfitserror(status);	
	
    fits_write_col(fptr,TDOUBLE,1,1,1,Nhyperparms,hypparms,&status);
	
	
    if (fits_close_file(fptr,&status))                
        printfitserror(status);           
	
}
