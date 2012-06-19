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
//  integrate_sed.c
//
//  Calculate the model flux for given filter.
//

#include "HBSGsep.h"

void integrate_sed(long length,double *sedval,double *fillam, \
				   double *filval,double *modval) {
	
    long ii;
    double h,w;
	
    h = *(fillam + 1) - *(fillam + 0);
	
    *modval=0.0;
    for (ii=0;ii < length-1;ii++) {
        if (ii == 0 || ii == length-2)
			w=3.0/8.0;
		else if (ii == 1 || ii == length-3)
			w=7.0/6.0;
		else if (ii == 2 || ii == length-4)
			w=23.0/24.0;
		else w=1.0;
		
		*modval = *modval + w * h * *(fillam+ii) * *(filval+ii) * \
				  *(sedval+ii);
    }
}
