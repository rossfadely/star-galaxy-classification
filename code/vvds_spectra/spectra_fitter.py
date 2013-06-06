import os
import commands
import ctypes
import numpy as np
import pyfits as pf

from ctypes import c_double,c_long,c_int,POINTER

class VVDSfitter(object):
    """
    Determine maximum likelihood classifications of VVDS spectra
    given template stars and galaxies
    """
    def __init__(self,wavelengths,spectradir,
                 uncertaintydir,starsdir,galsdir,
                 outname='VVDStemplatefits.dat',
                 dz=0.05,zmax=4.0):
        """
        Args are disk locations of directories containing only the
        appropriate VVDS FITS files, or spectral templates...
        except for wavelengths, which is the location of the
        wavelength file.
        """
        self.dz = dz
        self.zmax = zmax
        self.wavelengths = np.loadtxt(wavelengths)

        # Load data
        l = commands.getoutput('ls '+spectradir)
        self.spectra = self.get_spectra(l.split(),spectradir)
        l = commands.getoutput('ls '+uncertaintydir)
        self.uncertainty = self.get_spectra(l.split(),uncertaintydir)


        # Load templates
        l = commands.getoutput('ls '+starsdir)
        self.input_stars, self.Nstar = self.get_2col_ascii_data(l.split(),starsdir)
        l = commands.getoutput('ls '+galsdir)
        self.input_gals, self.Ngal = self.get_2col_ascii_data(l.split(),galsdir)
        
        self.stars = np.zeros((self.Nstar,self.wavelengths.shape[0]))
        for i in range(self.Nstar):
            self.stars[i,:] = self.input_stars[i][1,:]
        self.gals = np.zeros((self.Ngal,self.wavelengths.shape[0]))
        for i in range(self.Nstar):
            self.gals[i,:] = self.input_gals[i][1,:]
        # Redshifts galaxies
        #self.input_gals = self.redshift_gals()

        # Regrid templates
        #self.stars = self.regrid_templates(self.input_stars,self.Nstar)
        #self.gals = self.regrid_templates(self.input_gals,self.Ngal)

        # fit spectra
        self.star_scale,self.star_chi2,self.star_idx = self.fit_spectra('stars')
        self.gal_scale,self.gal_chi2,self.gal_idx = self.fit_spectra('gals')

        # write results
        self.write(outname)
        
    def get_2col_ascii_data(self,list,dir):
        """
        Return dictionary of wavelengths, flambdas for the list of star/gal
        templates
        """
        curdir = commands.getoutput('pwd')
        os.chdir(dir)
        templates = {}
        for i,l in enumerate(list):
            d = np.loadtxt(l)
            assert d.shape[1]==2
            templates[i] = d.T 
            if i%100==0: print 'Read %dth model' % i

        os.chdir(curdir)
        return templates, len(list)

    def get_spectra(self,list,dir):
        """
        Return numpy array of Ndata x Nwavelengths
        """
        curdir = commands.getoutput('pwd')
        os.chdir(dir)
        for i,l in enumerate(list):
            f = pf.open(l)
            d = f[0].data
            f.close()
            if i==0:
                spectra = np.zeros((len(list),d.shape[1]))
            spectra[i,:] = d[0,:] 
            if i%1000==0: print 'Read %dth spectrum' % i

        os.chdir(curdir)
        return spectra

    def redshift_gals(self):
        """
        Blow up galaxy grids in redshift space
        """
        newdict = {}
        zgrid = np.linspace(0.,self.zmax,self.zmax/self.dz+1)
        Nz = zgrid.shape[0]


        # again brutal foo
        for i in range(self.Ngal):
            for j in range(zgrid.shape[0]):
                d = self.input_gals[i].copy()
                d[0,:] *= (1.0 + zgrid[j])
                newdict[i*Nz+j] = d

        self.Ngal *= Nz
        return newdict

    def regrid_templates(self,template_dict,Ntemplates):
        """
        Regrid the templates to be the average in the bin
        """
        wVVDS = self.wavelengths


        # All VVDS spectra are on same grid
        Nwave = wVVDS.shape[0]
        hstep = (wVVDS[1] - wVVDS[0]) / 2.0
        regrided = np.zeros((Ntemplates,Nwave))

        wavep = wVVDS.ctypes.data_as(POINTER(c_double))
        mask = np.zeros(Nwave).astype('int32')
        mask_p = mask.ctypes.data_as(POINTER(c_int))

        interp = self._load_cubic_spline_interp_1d('./_cubic_spline_interp_1d.so',
                                                   'cubic_spline_interp_1d')

        # brutal, this is a bag of butts
        for i in range(Ntemplates):
            print 'regriding',i
            fnew = wVVDS * 0.0
            w = template_dict[i][0,:].astype('float64')
            f = template_dict[i][1,:].astype('float64')
            w_p = w.ctypes.data_as(POINTER(c_double))
            f_p = f.ctypes.data_as(POINTER(c_double))
            fnew_p = fnew.ctypes.data_as(POINTER(c_double))
            interp(w.shape[0],wVVDS.shape[0],w_p,f_p,wavep,fnew_p,mask_p)
            regrided[i,:] = fnew
            for j in range(Nwave):
                ind = (w>wVVDS[j]-hstep) & \
                    (w<=wVVDS[j]+hstep)
                if np.any(ind):
                    regrided[i,j] = np.sum(f[ind]) / float(f[ind].shape[0])
            if i%100==0: print 'Regrided the %dth template' % i
        return regrided


    def _load_cubic_spline_interp_1d(self,dll_path,function_name):
        """
        This reads in the compiled interpolation library
        """
        dll = ctypes.CDLL(dll_path,mode=ctypes.RTLD_GLOBAL)
        func = dll.cubic_spline_interp_1d
        func.argtypes = [c_long,c_long,POINTER(c_double),
                         POINTER(c_double),POINTER(c_double),
                         POINTER(c_double),POINTER(c_int)]
        return func

    def fit_spectra(self,type):
        """
        Fit template to spectra
        """
        if type=='stars':
            models = self.stars
        if type=='gals':
            models = self.gals

        Nspectra = self.spectra.shape[0]
        minchi2 = np.zeros(Nspectra)
        scales = np.zeros(Nspectra)
        template = np.zeros(Nspectra)

        # shoot me again
        for i in range(Nspectra):
            s = np.zeros(models.shape[0])
            c = np.zeros(models.shape[0])
            ind = self.uncertainty[i,:]>0.0
            iv = 1. / self.uncertainty[i,ind]**2.
            y  = self.spectra[i,ind]
            for j in range(models.shape[0]):
                m = models[j,ind]
                s[j] = np.dot(m,y * iv) / \
                    np.dot(m,m *iv)
                c[j] = np.sum((y-s[j]*m)**2.*iv)
            ind = np.argsort(c)[0]
            scales[i] = s[ind]
            minchi2[i] = c[ind]
            template[i] = ind
            if i%500==0: print 'Done Fitting the %dth spectrum' % i
        return scales, minchi2, template

    def write(self,outname):

        f = open(outname,'w')
        f.write('# starchi2 starscale starindex galchi2'+
                ' galscale galindex meanSN median SN\n')
        for i in range(self.spectra.shape[0]):
            ind = self.uncertainty[i,:]>0
            string = '%e %e %d %e %e %d %e %e\n' % (self.star_chi2[i],
                                                    self.star_scale[i],
                                                    self.star_idx[i],
                                                    self.gal_chi2[i],
                                                    self.gal_scale[i],
                                                    self.gal_idx[i],
                                                    np.mean(self.spectra[i,ind] / self.uncertainty[i,ind]),
                                                    np.median(self.spectra[i,ind] / self.uncertainty[i,ind]))
            if i%1000==0: print 'Wrote the %dth spectrum' % i
            f.write(string)
        f.close()
