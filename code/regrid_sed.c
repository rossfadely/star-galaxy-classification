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
//  regrid_sed.c
//
//  Regrid input SED profile to a common, finer wavelength grid that 
//  both SEDs and filters will share.
//

#include "HBSGsep.h"

void regrid_sed(double z,double *plam,double *pval,long finelength, \
				long sedlength,double *fineplam,double *finepval) {
	
    long ii;
    double x[sedlength],y[sedlength];
	
    for (ii=0;ii<sedlength;ii++) {
        x[ii]=*(plam+ii) * (1.0 + z);
	y[ii]=*(pval+ii);
    }
	
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline    = gsl_spline_alloc(gsl_interp_cspline, sedlength);
    
    gsl_spline_init(spline, x, y, sedlength);
	
    for (ii=0; ii < finelength; ii++) {
        if (*(fineplam+ii)/(1.0 + z) < *(plam+0)) {
	    *(finepval+ii) = 0.0;
	} else {
	    *(finepval+ii) = gsl_spline_eval (spline, *(fineplam+ii), acc);
	}
	if (*(finepval+ii) < 0.0)
	    *(finepval+ii) = 0.0;
    }
	
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}

