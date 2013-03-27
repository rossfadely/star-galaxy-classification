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
//  regrid.c
//
//  Regrid input SED/filter profile to a common, uniform wavelength grid 
//  that both SEDs and filters will share.
//


#include <stdio.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>


void regrid(double z,double wave_min,long oNstep,long nNstep,double step,
            double *owave,double *oval,double *nwave,double *nval) {
	
    long ii;
    double x[oNstep],y[oNstep];
	
    // redshift wavelengths
    for (ii=0;ii<oNstep;ii++) {
        x[ii] = owave[ii] * (1.0 + z);
        y[ii] = oval[ii];
    }
	
    // setup GSLs interpolation
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline    = gsl_spline_alloc(gsl_interp_cspline, oNstep);
    gsl_spline_init(spline, x, y, oNstep);

    // interpolate
    for (ii=0; ii < nNstep; ii++) {
        nwave[ii] = wave_min + (double)ii * step;
        if (nwave[ii] < owave[0] * (1.0 + z)) {
            nval[ii] = 0.0;
        } else {
            nval[ii] = gsl_spline_eval(spline, nwave[ii], acc);
        }
        if (nval[ii] < 0.0)
            nval[ii] = 0.0;
    }
	
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}

