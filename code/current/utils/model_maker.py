import subprocess
import numpy as np
import ctypes as ct

class ModelMaker(object):
    """
    Functionality to return model fluxes (and magnitudes), given a
    list of filters and SEDs.
    """
    def __init__(self,filter_dir=None,sed_dir=None,
                 z_grid=None,dll='./lib/_model_maker.so'):

        # get filters/seds
        self.sed_waves,self.sed_fluxes = self.get_list_and_read(sed_dir)
        self.filter_waves,self.filter_thrus = self.get_list_and_read(filter_dir)

        # calculate models
        self.make_models(dll,z_grid)

    def get_list_and_read(self,loc):
        """
        Use subprocess to get list of file names, then read
        the text files.  Return numpy arrays of wavelength and
        flux/throughput for seds/filters.
        """
        # optional manual input
        if loc==None:
            print 'Enter directory for SED/filter files:\n'
            loc = raw_input('>')

        # checking
        assert isinstance(loc,str), 'Directory entry should be str'
        if loc[-1]!='/':
            loc = loc+'/'

        # get list
        process = subprocess.Popen(['ls',loc], stdout=subprocess.PIPE)
        files, err = process.communicate()
        files = files.split()

        # create list
        wave = len(files)*[None]
        vals = len(files)*[None]
        for i in range(len(files)):
            tmp = np.loadtxt(loc+files[i])
            wave[i] = tmp[:,0] 
            vals[i] = tmp[:,1] 

        return wave,vals

    def make_models(self,dll,z_grid):
        """
        Produce model fluxes/mags
        """
        # initialize 
        self.Nseds = len(self.sed_waves)
        self.Nfilters = len(self.filter_waves)
        self.model_fluxes = np.zeros((Nseds,Nfilters)) 
        if z_grid==None:
            z_grid = np.zeros(1)

        # Ctypes foo
        maker = ct.CDLL(dll)
        sed_waves_pt = make_pointer_array(self.sed_waves)
        sed_fluxes_pt = make_pointer_array(self.sed_fluxes)
        filter_waves_pt = make_pointer_array(self.filter_waves)
        filter_thrus_pt = make_pointer_array(self.filter_thrus)
        model_fluxes_pt = make_pointer_array(self.model_fluxes)

    def make_pointer(self,arr):

    def make_pointer_array(self,arr):
        """
        Prep arrays for ctypes
        """
        ctarr = len(arr)*[None]
        for i in range(len(arr)):
            ctarr[i] = np.ctypeslib.as_ctypes(arr[i])

        return (ct.POINTER(ct.c_double)*len(arr))(*ctarr)


if __name__=='__main__':
    m = ModelMaker()
    print m.seds[0]
