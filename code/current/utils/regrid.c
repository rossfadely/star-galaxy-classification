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
//  Regrid input SED/filter profile to a common, finer wavelength grid 
//  that both SEDs and filters will share.
//


#include <stdio.h>
#include <gsl/gsl_spline.h>

void regrid(long, long, double, double*, double*, double*, double*);

void regrid(long clength, long flength, double z, double *cwave, \
            double *cflux, double *fwave, double *fflux) {
	
    long i;
    double x[clength],y[clength];
	
    // redshift wavelengths
    for (i=0;i<clength;i++) {
        x[i] = cwave[i] * (1.0 + z);
        y[i] = cflux[i];
    }
	
    // setup GSLs interpolation
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline    = gsl_spline_alloc(gsl_interp_cspline, clength);
    gsl_spline_init(spline, x, y, clength);

    // interpolate
    for (i=0; i < flength; i++) {
        if (fwave[i]/(1.0 + z) < cwave[0]) {
            fflux[i] = 0.0;
        } else {
            fflux[i] = gsl_spline_eval(spline, fwave[i], acc);
        }
        if (fflux[i] < 0.0)
            fflux[i] = 0.0;
    }
	
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}

