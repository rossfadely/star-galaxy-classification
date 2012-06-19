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
//  regrid_filter.c
//
//  Regrid input filter profile to a common, finer wavelength grid that 
//  both SEDs and filters will share.
//

#include "HBSGsep.h"

void regrid_filter(double *plam,double *pval,long slength,long length, \
				   double *fineplam,double *finepval) {
	
	long   ii;
	double min,max,step;
	double x[slength],y[slength];
	
	// original values
	for (ii=0;ii<slength;ii++) {
		x[ii] = *(plam+ii);
		y[ii] = *(pval+ii);
	}
	
	// initialize interpolation
	min  = x[0];
	max  = x[slength-1];
	step = (max-min)/(length-1);
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, slength);
	
	gsl_spline_init(spline, x, y, slength);
	
	// interpolate using spline, force throughput > 0
	for (ii=0; ii < (length); ii++) {
		*(fineplam+ii) = step * ii + min;
		*(finepval+ii) = gsl_spline_eval (spline, *(fineplam+ii), acc);
		
		if (*(finepval+ii) < 0.0)
			*(finepval+ii) = 0.0;
	}
	
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
}


