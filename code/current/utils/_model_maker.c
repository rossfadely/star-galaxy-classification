//
// This is part of HBSGsep, a hierarchical routine to classify astronomical 
// sources using photometric data.
//
// Please consult http://github.com/rossfadely/star-galaxy-classifiction
// for associated documentation.
//
//
// ==========================================================================
//
// model_maker.c
//
// Routine to produce model fluxes, given a list of SEDs and filters.  
// Also returned are the filter's normalization (zeropoint) in flux 
// density (F_lambda), for AB photometric system.
//
// ::filterlist:: - string specifying the path to filter list text file.
// ::sedlist:: - string specifying the path to sed list text file.
// ::Nz:: - Number of steps to ::zmax::, must be >=1.
// ::zmax:: - Maximum redshift for the models.
// ::norm:: - Pointer to array of filter normaliztions, must be same 
//            size array as number of filters. Filled out during runtime.
// ::models:: - Array of pointers to model fluxes.  Points to a MxN array 
//              where M is number of seds times Nz, and N is the number of 
//              filters.



#include <stdio.h>

void model_maker(char *filterlist,char *sedlist,long Nz,double zmax,
                 double *norm,double **models) {
	
    // Check Nz
    if (Nz<1) {
        printf("Nz must be greater than or equal to 1.\n\n");
        return;
    }
    
    // Declarations
    long ii,jj,kk;
    long *N,Nstep;
    long Nsed,Nfilter;
    long *filt_length;
    double *modval;
    double zgrid[Nz];
    double step,sed_min_step;
    double wave_min,wave_max;
    double *wave,*flux,*pnorm;
    double *sed_wave,*sed_flux;
    double *wave_new,*thru_new,*flux_new;
    double *filt_min_step;
    double **filt_lamb,**filt_thru;
    
    // Alloc
    N     = (long *)malloc(sizeof(long));
    pnorm = (double *)malloc(sizeof(double));
    
    // Create redshift grid
    if (Nz>1) {
        step = zmax / ((double)Nz - 1.0);
        for (ii=0; ii<Nz; ii++) {
            zgrid[ii] = step * (double)ii;
        }
    } else {
        zgrid[0] = 0.0;
    }

    
    // How many filters?
    get_num_files(filterlist, N);
    Nfilter = *N;
    if (Nfilter<2) {
        printf("Need more than 1 filter\n");
        return;
    }
    printf("\n\nFound %ld Filters\n",Nfilter);
    
    // Alloc
    filt_lamb = malloc(Nfilter*sizeof(double*));
    filt_thru = malloc(Nfilter*sizeof(double*));
    filt_length = (long *)malloc(Nfilter*sizeof(long));
    filt_min_step = (double *)malloc(Nfilter*sizeof(double));    
    
    // Read in the number of SEDs
    get_num_files(sedlist, N);
    Nsed = *N;
        
    // Read in filters, record smallest wavelength step, calc filter norm
    for (ii=0; ii<Nfilter; ii++) {
        get_filelength(ii,filterlist, N);
        filt_length[ii] = *N;
        filt_lamb[ii] = (double *)malloc(filt_length[ii] * sizeof(double));
        filt_thru[ii] = (double *)malloc(filt_length[ii] * sizeof(double));
        read_file(ii,filterlist,filt_lamb[ii],filt_thru[ii]);
        calc_normalization(filt_lamb[ii],filt_thru[ii], \
                           filt_length[ii],pnorm);
        norm[ii] = *pnorm;
        filt_min_step[ii] = filt_lamb[ii][1]-filt_lamb[ii][0];
        for (jj=1; jj<filt_length[ii]-1; jj++) {
            if ((filt_lamb[ii][jj+1]-filt_lamb[ii][jj])<filt_min_step[ii]) 
                filt_min_step[ii] = filt_lamb[ii][jj+1]-filt_lamb[ii][jj];
        }
    }

    // Read in SEDs, regrid, calculate fluxes
    for (ii=0; ii<Nsed; ii++) {
        // Read SEDs
        get_filelength(ii,sedlist, N);
        sed_wave = (double *)malloc(*N * sizeof(double));
        sed_flux = (double *)malloc(*N * sizeof(double));
        read_file(ii,sedlist,sed_wave,sed_flux);
        sed_min_step = sed_wave[1]-sed_wave[0];
        for (jj=1; jj<*N-1; jj++) {
            if ((sed_wave[jj+1]-sed_wave[jj])<sed_min_step) 
                sed_min_step = sed_wave[jj+1]-sed_wave[jj];
        }
        // Loop over redshift
        for (jj=0; jj<Nz; jj++) {
            // Loop over filter
            for (kk=0; kk<Nfilter; kk++) {
                wave_max = filt_lamb[kk][filt_length[kk]-1];
                wave_min = filt_lamb[kk][0];
                if (filt_min_step[kk]<sed_min_step) {
                    step = filt_min_step[kk];
                } else {
                    step = sed_min_step;
                }
                Nstep = (wave_max-wave_min) / (step) + 1;
                modval = (double *)malloc(sizeof(double));
                wave_new = (double *)malloc(Nstep * sizeof(double));
                thru_new = (double *)malloc(Nstep * sizeof(double));
                flux_new = (double *)malloc(Nstep * sizeof(double));
                regrid(0.0,wave_min,filt_length[kk],Nstep,step,filt_lamb[kk],
                       filt_thru[kk],wave_new,thru_new);
                regrid(zgrid[jj],wave_min,*N,Nstep,step,sed_wave,sed_flux,
                       wave_new,flux_new);
                integrate_sed(Nstep,flux_new,wave_new,thru_new,modval);
                models[ii*Nz+jj][kk] = *modval;
                free(modval);
                free(wave_new);
                free(thru_new);
                free(flux_new);
            }
        }
        free(sed_wave);
        free(sed_flux);
    }

    
    // Clean up
    for (ii=0; ii<Nfilter; ii++) {
        free(filt_lamb[ii]);
        free(filt_thru[ii]);
    }
    free(N);
    free(pnorm);
    free(filt_lamb);
    free(filt_thru);
    free(filt_min_step);
}


