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
//  get_num_files.c
//
//  Return number of files in an ascii list.
//

#include "HBSGsep.h"

void get_num_files(char *filelocations,long *N) {
		
	FILE *fp;
	char filepath[3000];
	
	*N = 0;
	fp = fopen(filelocations,"r");
	while (fgets(filepath,FILEPATH_LENGTH,fp) != NULL) {
		*N += 1;
	}	
	fclose(fp);
}

