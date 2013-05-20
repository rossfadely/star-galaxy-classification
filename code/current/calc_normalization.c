//
// This is part of HBSGsep, a hierarchical routine to classify stars and 
// galaxies using photometric data.
//
// Please consult http://github.com/rossfadely/star-galaxy-classifiction
// for associated documentation.
//
//
// ==========================================================================
//
//  calc_normalization.c
//
//  Calculate the specified filter's normalization (zeropoint) in flux 
//  density (F_lambda), for AB photometric system.
//

#include <stdio.h>

void calc_normalization(double *lam,double *val,long length, \
                        double *normval) {

    double abfnu = 3.631e-20, cinang = 3.0e18;
    double tmpint, h, w;
    long ii;

    tmpint=0;
    h = *(lam+1) - *(lam +0);
	
    for (ii=0;ii < length-1;ii++) {
        if (ii == 0 || ii == length-2)
            w=3.0/8.0;
        else if (ii == 1 || ii == length-3)
            w=7.0/6.0;
        else if (ii == 2 || ii == length-4)
            w=23.0/24.0;
        else w=1.0;
		
        tmpint += w * h / lam[ii] * val[ii] * abfnu * cinang;
		
    }
    *normval = tmpint;
}


