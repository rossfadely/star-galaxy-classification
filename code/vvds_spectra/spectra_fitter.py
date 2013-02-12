import commands
import ctypes
import numpy as np
import pyfits as pf

class VVDSfitter(object):
    """
    Determine maximum likelihood classifications of VVDS spectra
    given template stars and galaxies
    """
    def __init__(self,spectra,uncertainty,stars,gals,dz=0.05,zmax=4.0):
        """
        Args are disk locations of directories containing only the
        appropriate VVDS FITS files, or spectral templates
        """
        
        # Get lists of files
        l = commands.getoutput('ls '+spectra)
        self.spectra_files = l.split()
        l = commands.getoutput('ls '+uncertainty)
        self.uncertainty_files = l.split()
        l = commands.getoutput('ls '+stars)
        self.stars_files = l.split()
        l = commands.getoutput('ls '+gals)
        self.gals_files = l.split()

    def get_2col_ascii_data(self,list):
        """
        Return dictionary of wavelengths, flambdas for the list of star/gal
        templates
        """
        flambdas = {}
        wavelengths = {}
        for i,l in enumerate(list):
            d = np.loadtxt(l)
            assert d.shape[0]==2
            wavelengths[i] = d[0,:] 
            flambdas[i] = d[1,:] 

        return wavelengths, flambdas
