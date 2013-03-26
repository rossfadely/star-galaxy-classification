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
//  calc_AB_normalization.c
//
//  Calculate the specified filter's normalization (zeropoint) in flux 
//  density (F_lambda), for AB photometric system.
//

#include <stdio.h>

void calc_AB_normalization(double *wave,double *thru,double length, \
                           double *normval) {

    double abfnu = 3.631e-20, cinang = 3.0e18;
    double tmpint, h, w;
    long i;

    tmpint=0;
    h = wave[1] - wave[0];
	
    for (i=0;i < length-1;i++) {
        if (i == 0 || i == length-2)
            w=3.0/8.0;
        else if (i == 1 || i == length-3)
            w=7.0/6.0;
        else if (i == 2 || i == length-4)
            w=23.0/24.0;
        else w=1.0;
		
        tmpint += w * h / wave[i] * flux[i] * abfnu * cinang;
		
    }
    *normval = tmpint;
}


