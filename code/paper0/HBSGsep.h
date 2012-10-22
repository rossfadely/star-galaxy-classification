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
// HBSGsep.h
//
// This is the include file needed to run the program - it declares 
// needed includes as well as many global variables used in different 
// routines.  GSL is needed to run the program, the library is available 
// at http://www.gnu.org/software/gsl/ . Also needed are the CFITSIO 
// routines, available at http://heasarc.gsfc.nasa.gov/fitsio/ .
//





// Standard includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

// GSL includes
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>

// CFITSIO include
#include "fitsio.h"

// change if length of any file path > 3000 characters
#define FILEPATH_LENGTH 3000


//storage variables
long Nfilter,Ndata;
long Nstartemplate,Ngaltemplate;
long Nstarhyperparms,Ngalhyperparms;
long Nhyperparms;
long count,count_tot;
long doneflag;
long *filter_lgth_fine; 
long *datatype;
long *zinds;
long **starind,**galind;
double zstep;
double lnPtot_start;
double old_lnPtot;
double size;
double *norm;
double *lnlikeratio;
double *star_coeff_mean,*star_coeff_var;
double *gal_coeff_mean,*gal_coeff_var;
double *star_minchi,*gal_minchi;
double **filter_lamb_fine, **filter_thru_fine;
double **dataflux, **datafluxerr;
double **modelflux_stars,**modelflux_gals;
double **galinfo,**starinfo;
double ***starsparse,***galsparse;




// Routines in order they are called
void filter_calcs(void);
void get_num_files(char*,long*);
void get_filelength(long,char*,long*);
void read_file(long,char*,double*,double*);
void regrid_filter(double*,double*,long,long,double*,double*);
void calc_normalization(double*,double*,double,double*);
void data_calcs(void);
void count_data(char*,long*);
void calc_datavals(long,char*,double,double*,double*);
void printfitserror(int status);
void star_model_calcs(void);
void regrid_sed(double,double*,double*,long,long,double*,double*);
void integrate_sed(long,double*,double*,double*,double*);
void star_model_mags(void);
void write_modelmags(char*,long,double**);
void gal_model_calcs(void);
void gal_model_mags(void);
void fit_calcs(void);
void coeff_calcs(void);
void fit_star_template(long,double*,double*,double*);
void coeff_meanvar(int,long,double*,double*,double*);
void fit_gal_template(long,long,double*,double*,double*);
void sparse_calcs(void);
void calc_P_F_kS(long);
void calc_P_F_kG(long);
void optimize(void);
void read_hyppars(char*,double*);
void calc_P_F_S(long,double*,double*);
void calc_P_F_G(long,double*,double*);
void write_lnprob(long,double);
void write_hyppars(char*,double*);
void write_lnlikeratio(char*,double*);
void write_minchi2(char*,char*,long);
double loglikelihood(const gsl_vector*,double*); 


