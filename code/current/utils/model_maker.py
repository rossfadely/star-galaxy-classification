import subprocess
import numpy as np

class ModelMaker(object):
    """
    Functionality to return model fluxes (and magnitudes), given a
    list of filters and SEDs.
    """
    def __init__(self,filter_dir=None,sed_dir=None,
                 dll_dir='./lib/',z_grid=None,
                 calc_norm=True):

        self.seds = self.get_list_and_read(sed_dir)
        self.filters = self.get_list_and_read(filter_dir)

    def get_list_and_read(self,loc):
        """
        Use subprocess to get list of file names, then read
        the text files.  Return list of 2D numpy arrays for
        filters/seds.
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
        tmp = len(files)*[None]
        for i in range(len(files)):
            tmp[i] = np.loadtxt(loc+files[i])

        return tmp






if __name__=='__main__':
    m = ModelMaker()
    print m.seds[0]
