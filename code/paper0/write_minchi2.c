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
// write_minchi2.c
//
// Write out minimum chi2 values. Type = 0 prints galaxy, else prints 
// star.
//
//

#include "HBSGsep.h"

void write_minchi2(char *outfile,char *columnname,long type) {
	
	
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
    strcpy(*(ttype+0),columnname);
    strcpy(*(tform+0),"1D");
    strcpy(*(tunit+0),"      ");
	
    remove(outfile);   

    if (fits_create_file(&fptr,outfile,&status)) 
        printfitserror(status);           
	
    if (fits_create_img(fptr,bitpix,naxis,naxes,&status))
        printfitserror(status);   
	
    if (fits_movabs_hdu(fptr,1,&hdutype,&status))
        printfitserror(status);
	
    if (fits_create_tbl(fptr,BINARY_TBL,Ndata,tfields,ttype,tform,tunit, \
                        extname,&status))
        printfitserror(status);	
	
    if (type == 0) {
        fits_write_col(fptr,TDOUBLE,1,1,1,Ndata,gal_minchi,&status);
    } else {
        fits_write_col(fptr,TDOUBLE,1,1,1,Ndata,star_minchi,&status);	
    }
	
		
    if (fits_close_file(fptr,&status))                
        printfitserror(status);           
	

}
