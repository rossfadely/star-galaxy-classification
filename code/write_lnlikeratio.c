/*
 *  write_lnlikeratio.c
 *  
 *	Write out ln(likelihood) ratios at end.
 *
 *  Created by Ross Fadely on 6/29/11.
 *  Copyright 2011 Haverford College. All rights reserved.
 *
 */

#include "HBSGsep.h"

void write_lnlikeratio(char *lnlikeratiooutfile,double *lnlikeratio) {
	
	
	// Initialize fits file
	int bitpix = USHORT_IMG; 
	long naxis = 1;
	long naxes[1] = {1};
	int status=0,tfields=1,hdutype;
    char *ttype[1],*tform[1],*tunit[1];          
    char extname[] = " ";           
	fitsfile *fptr;       
		
	ttype[0] = (char *)malloc(50 * sizeof(char));
	tform[0] = (char *)malloc(50 * sizeof(char));
	tunit[0] = (char *)malloc(50 * sizeof(char));
	strcpy(*(ttype+0),"ln_P_S_over_P_G");
	strcpy(*(tform+0),"1D");
	strcpy(*(tunit+0),"      ");
	
	remove(lnlikeratiooutfile);   

	if (fits_create_file(&fptr,lnlikeratiooutfile,&status)) 
		printfitserror(status);           
	
    if (fits_create_img(fptr,bitpix,naxis,naxes,&status))
		printfitserror(status);   
	
	if (fits_movabs_hdu(fptr,1,&hdutype,&status))
		printfitserror(status);
	
	if (fits_create_tbl(fptr,BINARY_TBL,Ndata,tfields,ttype,tform,tunit,extname,&status))
		printfitserror(status);	
	
	fits_write_col(fptr,TDOUBLE,1,1,1,Ndata,lnlikeratio,&status);
		
    if (fits_close_file(fptr,&status))                
		printfitserror(status);           
	

}
