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
// HBSGsep.c
//
// The main routine from which functions are called. The first portion 
// of the document lists variables the user *must* specify, and some 
// which can be modified to alter performance.  Please see documentation
// for details about user specified variables.
//

#include "HBSGsep.h"


/***********************************************************************
 User specified variables
 **********************************************************************/


// MUST CHANGE THE FOLLOWING
// *************************


// Input Files
char filtersinput[]  = "YOUR_PATH";
char galssedinput[]  = "YOUR_PATH";
char starssedinput[] = "YOUR_PATH";
 
// Data File 
char datafile[] = "YOUR_PATH";
 
// Output Files
char hypoutfile[]         = "YOUR_PATH";
char lnlikeratiooutfile[] = "YOUR_PATH";
char lnPtotfile[]         = "YOUR_PATH";

// Analysis variables
double zmin       = 0.0;
double zmax       = 4.0;
long Nz           = 51;

double noisefudge = 0.06;

long Niter        = 50000;													
float tol         = 1.0e-1;													


// CAN CHANGE THE FOLLOWING
// ************************


// Input Files
char hypinfile[]     = "YOUR_PATH";

// Output Files
char starmodmagsfile[]    = "YOUR_PATH";
char galmodmagsfile[]     = "YOUR_PATH";
char starchi2file[]       = "YOUR_PATH";
char galchi2file[]        = "YOUR_PATH";


// Analysis variables
double fluxunitfactor = 5.0e13;										
double probfrac       = 1e-10;												
double P_floor        = 1.0e-300;										    
double weight_fact    = 1.e-1;
long Ncstep           = 7;													   
long writeiter        = 1000;


// Analysis flags
long usehypin        = 0;		     				
long calc_model_mags = 0;
long writeminchi2    = 0;




/***********************************************************************
 Do HBSGsep
 **********************************************************************/
int main(int argc, char **argv) {
	
	long ii,jj;
	
	// Setup timing
	time_t  t0; 
	t0 = time(NULL);
	
	// Do Filter Calcs	
	filter_calcs();
	
	// Read Data, and do calculate fluxes, etc.
	data_calcs();
	
	// Calculate models for stars	
	star_model_calcs();
	
	// If desired, calculate star model mags
	if (calc_model_mags==1) {
		star_model_mags();
	}
	
	// Calculate models for gals	
	gal_model_calcs();
	
	// If desired, calculate gal model mags
	if (calc_model_mags==1) {
		gal_model_mags();
	}
	
	// Fit the templates to calculate the 
	// coefficient priors and make sparse 
	// arrays for optimizer.	
	fit_calcs();
	
	
	
	// Assign global number of hyperparameters
	Nstarhyperparms = Nstartemplate;
	Ngalhyperparms  = Ngaltemplate;
	Nhyperparms     = Nstarhyperparms+Ngalhyperparms+2;
	
	
	// Time to optimize.  Its useful to call 
	// this more than once, hence we read in 
	// the values of the hyperparameters from 
	// the preceding runs.
	if (Niter>0) {
		
		count_tot = 0;
		
	    optimize();
	    usehypin = 1;
	    optimize();
	    
		// Write the log likelihood ratios to file
	    write_lnlikeratio(lnlikeratiooutfile,lnlikeratio);
	}
	
	
	
	//  Write minimum chi2 values, if desired
	if (writeminchi2 == 1) {
		write_minchi2(starchi2file,"min_star_chi2",1);
		write_minchi2(galchi2file,"min_gal_chi2",1);
	}	
	
	
	// Cleanup allocated memory
	for (ii=0; ii<Nfilter; ii++) {
		free(filter_lamb_fine[ii]);
		free(filter_thru_fine[ii]);	
		free(dataflux[ii]);
		free(datafluxerr[ii]);
		free(modelflux_stars[ii]);
		free(modelflux_gals[ii]);
	}
	free(dataflux);
	free(datafluxerr);
	free(modelflux_stars);
	free(modelflux_gals);
	free(norm);
	free(filter_lamb_fine);
	free(filter_thru_fine);	
	free(filter_lgth_fine);	
	for (jj=0; jj<3; jj++) {
		for (ii=0; ii<Ndata; ii++) {
			free(starsparse[jj][ii]);
		}
		free(starsparse[jj]);
	}
	free(starsparse);
	for (jj=0; jj<3; jj++) {
		for (ii=0; ii<Ndata; ii++) {
			free(galsparse[jj][ii]);
		}
		free(galsparse[jj]);
	}
	free(galsparse);
	free(star_coeff_mean);
	free(star_coeff_var);
	free(star_minchi);
	free(gal_minchi);

	
	// Print elapsed time	
	printf("\n\nDONE!!!!\n\n");
	printf ("Elapsed time (s): %ld\n", (long) (time(NULL) - t0));
	
	return 0;
}









