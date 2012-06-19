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
// write_lnprob.c
//
// Write the current value of the loglikelihood to a simple text file 
// for monitoring purposes.
//
//

#include "HBSGsep.h"

void write_lnprob(long ii,double lnPtot) {
	
	extern char lnPtotfile[FILEPATH_LENGTH];

	FILE *fp;
	
	if (ii==0) {
		fp = fopen(lnPtotfile,"w");
		fprintf(fp,"\n");
		fclose(fp);
	}
	
	fp = fopen(lnPtotfile,"a+");  	
	fprintf(fp,"%ld %g\n",ii,lnPtot);
	fclose(fp);
}
