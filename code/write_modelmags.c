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
//  write_modelmags.c
//
//  Write model magnitudes to FITS file.  This is hard coded to ugriz.
//

#include "HBSGsep.h"

void write_modelmags(char *outfile,long Nmodels,double **modelmags) {
	
    // Initialize fits file
    int bitpix   =  USHORT_IMG; 
    long naxis    =   1;
    long naxes[1] = {1};
    int status=0,tfields=5,hdutype,kk;
    char *ttype[5],*tform[5],*tunit[5];          
    char extname[] = " ";           
    fitsfile *fptr;       
	
	
    ttype[0]=(char *)malloc(10 * sizeof(char));
    tform[0]=(char *)malloc(10 * sizeof(char));
    tunit[0]=(char *)malloc(10 * sizeof(char));
    strcpy(*(ttype+0),"u");
    strcpy(*(tform+0),"1D");
    strcpy(*(tunit+0),"      ");
    ttype[1]=(char *)malloc(10 * sizeof(char));
    tform[1]=(char *)malloc(10 * sizeof(char));
    tunit[1]=(char *)malloc(10 * sizeof(char));
    strcpy(*(ttype+1),"g");
    strcpy(*(tform+1),"1D");
    strcpy(*(tunit+1),"      ");
    ttype[2]=(char *)malloc(10 * sizeof(char));
    tform[2]=(char *)malloc(10 * sizeof(char));
    tunit[2]=(char *)malloc(10 * sizeof(char));
    strcpy(*(ttype+2),"r");
    strcpy(*(tform+2),"1D");
    strcpy(*(tunit+2),"      ");
    ttype[3]=(char *)malloc(10 * sizeof(char));
    tform[3]=(char *)malloc(10 * sizeof(char));
    tunit[3]=(char *)malloc(10 * sizeof(char));
    strcpy(*(ttype+3),"i");
    strcpy(*(tform+3),"1D");
    strcpy(*(tunit+3),"      ");
    ttype[4]=(char *)malloc(10 * sizeof(char));
    tform[4]=(char *)malloc(10 * sizeof(char));
    tunit[4]=(char *)malloc(10 * sizeof(char));
    strcpy(*(ttype+4),"z");
    strcpy(*(tform+4),"1D");
    strcpy(*(tunit+4),"      ");
	
    remove(outfile);   
	
    if (fits_create_file(&fptr,outfile,&status)) 
        printfitserror(status);	
	
    if (fits_create_img(fptr,bitpix,naxis,naxes,&status))
        printfitserror(status);	
	
    if (fits_movabs_hdu(fptr,1,&hdutype,&status))
        printfitserror(status);	
	
    if (fits_create_tbl(fptr,BINARY_TBL,Nmodels,tfields,ttype,tform,tunit, \
						extname,&status))
        printfitserror(status);	
	
    for (kk=0; kk<Nfilter; kk++) {
        fits_write_col(fptr,TDOUBLE,kk+1,1,1,Nmodels,modelmags[kk],&status);
    }
		
    if (fits_close_file(fptr,&status))                
        printfitserror(status);	
	
}
