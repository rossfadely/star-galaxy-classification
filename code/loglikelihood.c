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
// loglikelihood.c
//
// Return negative loglikelihood of entire data, given the current set 
// of hyperparameters.
//

#include "HBSGsep.h"

double loglikelihood(const gsl_vector *v, double *junk) {

    //external variables
	extern long writeiter;
    extern double old_lnPtot;
    extern char hypoutfile[FILEPATH_LENGTH];
	
    //internal variables
    double *hypparms_stars,*hypparms_gals,*hypparms_weights,hyptot;
    double *hypparms;
    double *P_F_S,*P_F_G;
    double lnPtot=0,lnPtot_G=0,lnPtot_S=0;
    long ii,ind;
	
    P_F_S = (double *)calloc(Ndata,sizeof(double));
    P_F_G = (double *)calloc(Ndata,sizeof(double));

    if (count_tot==0) {
		lnlikeratio = (double *)calloc(Ndata,sizeof(double));
    }
	
    hypparms_stars   = (double *)malloc(Nstarhyperparms * sizeof(double));
    hypparms_gals    = (double *)malloc(Ngalhyperparms * sizeof(double));
    hypparms_weights = (double *)malloc(2 * sizeof(double));
    hypparms         = (double *)malloc(Nhyperparms * sizeof(double));
	
    // Read in current value of hyperparms from optimizer
    // star hyperparms
	hyptot=0.0;
    for (ii=0; ii<Nstartemplate; ii++) {
        hypparms_stars[ii] = pow(10.0,gsl_vector_get(v,ii));
		hypparms_stars[ii] = sqrt(hypparms_stars[ii] * hypparms_stars[ii]);
		hyptot += hypparms_stars[ii];
    }
    for (ii=0; ii<Nstartemplate; ii++) {
        hypparms_stars[ii] /= hyptot;
		hypparms[ii] = hypparms_stars[ii];
    }
	
	// galaxy hyperparms
    hyptot=0.0;
    for (ii=Nstarhyperparms; ii<Nstarhyperparms+Ngalhyperparms; ii++) {
		ind = ii - Nstarhyperparms;
        hypparms_gals[ind] = pow(10.,gsl_vector_get(v,ii));
		hypparms_gals[ind] = sqrt(hypparms_gals[ind] * hypparms_gals[ind]);
		hyptot += hypparms_gals[ind];
    }
    for (ii=Nstarhyperparms; ii<Nstarhyperparms+Ngalhyperparms; ii++) {
		ind = ii - Nstarhyperparms;
        hypparms_gals[ind] /= hyptot;
		hypparms[ii] = hypparms_gals[ind];
    }
	
	// Weight hyperparms
    hyptot=0.0;
    for (ii=Nhyperparms-2; ii<Nhyperparms; ii++) {
		ind = ii - Nstarhyperparms - Ngalhyperparms;
        hypparms_weights[ind] = gsl_vector_get(v,ii);
		hypparms_weights[ind] = sqrt(hypparms_weights[ind] * \
									 hypparms_weights[ind]);
		hyptot += hypparms_weights[ind];
    }
    for (ii=Nhyperparms-2; ii<Nhyperparms; ii++) {
		ind = ii - Nstarhyperparms - Ngalhyperparms;
        hypparms_weights[ind] /= hyptot;
		hypparms[ii] = hypparms_weights[ind];
    }


    // Calculate P_F_S
    for (ii=0; ii<Ndata; ii++) {
		
        calc_P_F_S(ii,hypparms_stars,P_F_S);
		calc_P_F_G(ii,hypparms_gals,P_F_G);
		
		if (P_F_G[ii]==0)
			P_F_G[ii]=1.e-300;
		if (P_F_S[ii]==0)
			P_F_S[ii]=1.e-300;
		
		lnPtot_G += log(P_F_G[ii]);
		lnPtot_S += log(P_F_S[ii]);
		lnPtot   += log(P_F_S[ii] * hypparms_weights[0] + \
						P_F_G[ii] * hypparms_weights[1]);
				
		if (doneflag!=0) {
			lnlikeratio[ii] = log(P_F_S[ii] * hypparms_weights[0]) - \
							  log(P_F_G[ii] * hypparms_weights[1]);
		}
    }
	
    // Record starting point for the record
    if (count_tot==0) {
        lnPtot_start = lnPtot;
		old_lnPtot   = lnPtot;
    }
	
    // If the loglikelihood improves, write it out
    if (lnPtot > old_lnPtot) {
        write_hyppars(hypoutfile,hypparms);
		old_lnPtot   = lnPtot;
    }

	// Print status
	if (count_tot % writeiter == 0) {
		printf("\n\n\nTotal Count %ld\n",count_tot);		
		printf("Current Count %ld\n",count);		
		printf("Initial loglike %g\n",lnPtot_start);		
		printf("Current loglike %g\n",old_lnPtot);		
		printf("Simplex size %g\n",size);		
		fflush(stdout);
		sleep(1);
	}

    count     += 1;
    count_tot += 1;
	
    // Make it positive for minimizer
    lnPtot=sqrt(lnPtot*lnPtot);

    free(hypparms_stars);
    free(hypparms_gals);
    free(hypparms_weights);
    free(P_F_S);
    free(P_F_G);

    return lnPtot;
}
