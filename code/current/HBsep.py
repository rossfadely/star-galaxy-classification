
import numpy as np
import matplotlib.pyplot as pl

class HBsep(object):

    def __init__(self, data = None,
                 magcols = None,
                 errcols = None,
                 filterlist = None,
                 Nclasses = None):

        self.data = data
        self.read_data()
        self.filterlist = filterlist
        self.Nclasses   = Nclasses

    def read_data(self):
        """
        Read in data file or array, do some moderate
        error checking.
        """
        if self.data==None:
            print '\nSpecification for data not given,'
            print 'please give a file location or 2D numpy array'
            assert False
        elif isinstance(self.data,np.ndarray):
            self.shape_check()
            self.process_data()
        elif isinstance(self.data,str):
            try:
                # hdu number might be an issue
                f = pf.open(self.data)
                self.data = f[1].data
                f.close()
                self.shape_check()
                self.process_data()
            except:
                print '\nData file not read.'
                assert False

    def shape_check(self):
        """
        Simple shape checks
        """
        assert len(self.data.shape)==2, \
            '\nData must be a 2D numpy array'
        assert self.data.shape[0]>self.data.shape[1], \ 
            '\nNdata must be greater than Nfilter'
        assert np.mod(self.data.shape[1],2)==0, '\nNcolumn is not even,' + \
            'numpy array must have columns of mags then errs.'


    def process_data(self)
        """
        Turn data into mags, mag err, flux, flux err
        """
        self.Nfilter = self.data.shape[1] / 2
        if keepmags:
            self.mags = self.data[:,:Nfilter]
            self.merr = self.data[:,Nfilter:]

        # Treat missing data
