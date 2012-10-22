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
//  get_filelength.c
//
//  From list of file locations, read in the specified file and return 
//  its length.
//

#include "HBSGsep.h"

void get_filelength(long filenum,char *filelocations,long *N) {
	
    FILE   *fp,*fp2;
    char   filepath[3000];
    short  mod;
    long   ii = 0,jj = 0,kk = 0;
    double ftmp;
	
    fp = fopen(filelocations,"r");
    while (fgets(filepath,FILEPATH_LENGTH,fp) != NULL) {
        if (ii == filenum) {
            char *pchr = strchr(filepath, '\n');
            if(pchr != NULL) *pchr='\0';
            fp2 = fopen(filepath, "r");
            while (fscanf (fp2, "%lf", &ftmp) == 1) {
                mod = jj % 2;
                if (mod == 1) {
                    kk++;
                }
                jj++;
            }
            *N = kk;			
            fclose(fp2);
        }
        ii += 1;
    }	
    fclose(fp);
}
