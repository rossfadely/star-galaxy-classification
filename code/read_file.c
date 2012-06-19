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
//  read_file.c
//
//  Read in contents of two column ascii data.
//

#include "HBSGsep.h"

void read_file(long filenum,char *filelocations,double *plam,double *pval) {
	
    FILE   *fp,*fp2;
    char   filepath[FILEPATH_LENGTH];
    long   ii = 0,jj = 0,kk = 0;
    short  mod;
    double ftmp;
	
	
    fp = fopen(filelocations,"r");
    while (fgets(filepath,FILEPATH_LENGTH,fp) != NULL) {
        if (ii == filenum) {
            char *pchr = strchr(filepath, '\n');
            if(pchr != NULL) *pchr='\0';
            fp2 = fopen(filepath, "r");
            while (fscanf (fp2, "%lf", &ftmp) == 1) {
                mod = jj % 2;
                if (mod==0) {
                    *(plam + kk)=ftmp;
                } else {
                    *(pval + kk)=ftmp;
                    kk++;
                }
                jj++;
            }			
            fclose(fp2);
        }
        ii++;
    }	
    fclose(fp);	
}


