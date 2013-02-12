import os
import commands
import ctypes
import numpy as np
import pyfits as pf

class VVDSfitter(object):
    """
    Determine maximum likelihood classifications of VVDS spectra
    given template stars and galaxies
    """
    def __init__(self,wavelengths,spectradir,
                 uncertaintydir,starsdir,galsdir,dz=0.05,
                 zmax=4.0):
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

        # Calc S/N
        self.SNspectra = np.abs(self.spectra) / self.uncertainty

        # Load templates
        l = commands.getoutput('ls '+starsdir)
        self.input_stars, self.Nstar = self.get_2col_ascii_data(l.split(),starsdir)
        l = commands.getoutput('ls '+galsdir)
        self.input_gals, self.Ngal = self.get_2col_ascii_data(l.split(),galsdir)

        # Redshifts galaxies
        self.input_gals = self.redshift_gals()

        # Regrid templates
        self.stars = self.regrid_templates(self.input_stars,self.Nstar)

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
        return newdict



    def regrid_templates(self,template_dict,Ntemplates):
        """
        Regrid the templates to be the average in the bin
        """

        # All VVDS spectra are on same grid
        Nwave = self.wavelengths.shape[0]
        hstep = (self.wavelengths[1] - self.wavelengths[0]) / 2.0
        regrided = np.zeros((Ntemplates,Nwave))
        # NOT WORKING
        # brutal, this is a bag of butts
        for i in range(Ntemplates):
            w = template_dict[i][0,:]
            f = template_dict[i][1,:]
            for j in range(Nwave):
                ind = (w>self.wavelengths[j]-hstep) & \
                    (w<=self.wavelengths[j]-hstep)
                if np.any(ind):
                    regrided[i,j] = np.sum(f[ind]) / float(f[ind].shape[0])

        return regrided



        
        
