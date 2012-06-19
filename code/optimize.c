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
// optimize.c
//
// Set up optimizer, run it, and then report how it went.  The optimizer 
// calcs the function loglikelihood.c
//

#include "HBSGsep.h"

void optimize(void) {
	
    //External variables
    extern double lnPtot_start,weight_fact;
    extern double old_lnPtot;
    extern long Niter,writeiter;
    extern float tol;
    extern long usehypin;
    extern char hypinfile[FILEPATH_LENGTH],hypoutfile[FILEPATH_LENGTH];

    //Internal variables
    size_t iter = 0;
    int status,ii;
    double *hypparms,hyptot;
	
    hypparms=(double *)malloc((Nhyperparms) * sizeof(double));
    count    = 0;
    doneflag = 0;
		
    // GSL optimizer setup
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s =  NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function opt_func;
	
    // Initialize to 1/Ntemplate or read in existing file
    x = gsl_vector_alloc(Nhyperparms);
    if (usehypin!=0) {

        read_hyppars(hypinfile,hypparms);
		for (ii=0; ii<Nhyperparms-2; ii++) {
			gsl_vector_set(x,ii,log10(hypparms[ii]));
		}
		for (ii=Nhyperparms-2; ii<Nhyperparms; ii++) {
			gsl_vector_set(x,ii,hypparms[ii]);
		}
		
    } else {

		for (ii=0; ii<Nstarhyperparms; ii++) {
			gsl_vector_set(x,ii,log10(1.0/((double)Nstartemplate)));
		}
		for (ii=Nstarhyperparms; ii<Ngalhyperparms+Nstarhyperparms; ii++) {
			gsl_vector_set(x,ii,log10(1.0/((double)Ngaltemplate)));
		}
		for (ii=Nhyperparms-2; ii<Nhyperparms; ii++) {
			gsl_vector_set(x,ii,0.5);
		}	
		
    }
	
    // Initialize step sizes
    ss = gsl_vector_alloc(Nhyperparms);

    for (ii=0; ii<Nstarhyperparms; ii++) {
        gsl_vector_set(ss,ii,5.0);
    }

    for (ii=Nstarhyperparms; ii<Ngalhyperparms+Nstarhyperparms; ii++) {
        gsl_vector_set(ss,ii,5.0);
    }
    for (ii=Nhyperparms-2; ii<Nhyperparms; ii++) {
        gsl_vector_set(ss,ii,gsl_vector_get(x,ii)*weight_fact);
    }	
			
    // Setup function
    opt_func.n      = Nhyperparms;
    opt_func.f      = loglikelihood;
    opt_func.params = hypparms;   
	
    // Set optimizer loose
    s = gsl_multimin_fminimizer_alloc(T,Nhyperparms);
    gsl_multimin_fminimizer_set(s,&opt_func,x,ss);
    do {
		
		iter++;
		status=gsl_multimin_fminimizer_iterate(s);
		
		if (status)
			break;
		
		size   = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size,tol);
		if (status==GSL_SUCCESS) {
			doneflag = 1;
			status   = gsl_multimin_fminimizer_iterate(s);
		}
		
		if (count_tot % writeiter == 0) {
			write_lnprob((iter-1),s->fval);
		}
		
    } while (status==GSL_CONTINUE && iter < Niter);
    
    doneflag = 1;
    status   = gsl_multimin_fminimizer_iterate(s);

	
    // Its done, renormalize hyperparms that are spit out	
    hyptot=0.0;
    for (ii=0; ii<Nstarhyperparms; ii++) {
        hypparms[ii] = pow(10.,gsl_vector_get(s->x,ii));
		hypparms[ii] = sqrt(hypparms[ii] * hypparms[ii]);
		
		hyptot += hypparms[ii];
    }
    for (ii=0; ii<Nstarhyperparms; ii++) {
        hypparms[ii] /= hyptot;
    }
    hyptot=0.0;
    for (ii=Nstarhyperparms; ii<Ngalhyperparms+Nstarhyperparms; ii++) {
        hypparms[ii] = pow(10.,gsl_vector_get(s->x,ii));
		hypparms[ii] = sqrt(hypparms[ii] * hypparms[ii]);
		
		hyptot += hypparms[ii];
    }
    printf("\n");
    for (ii=Nstarhyperparms; ii<Ngalhyperparms+Nstarhyperparms; ii++) {
        hypparms[ii] /= hyptot;
    }
    hyptot=0.0;
    for (ii=Nhyperparms-2; ii<Nhyperparms; ii++) {
        hypparms[ii] = gsl_vector_get(s->x,ii);
		hypparms[ii] = sqrt(hypparms[ii] * hypparms[ii]);
		
		hyptot += hypparms[ii];
    }
    printf("\n");
    for (ii=Nhyperparms-2; ii<Nhyperparms; ii++) {
        hypparms[ii] /= hyptot;
    }
	
	
    //Show how we did
    printf("\n\nSummary:\n");
    printf("Initial loglikelihood: %g \n",lnPtot_start);
    printf("  Final loglikelihood: %g \n",old_lnPtot);
		
    //write the hyperparms to file
    write_hyppars(hypoutfile,hypparms);

    //cleanup
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
    
}
