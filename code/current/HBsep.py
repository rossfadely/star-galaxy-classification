
import numpy as np
import matplotlib.pyplot as pl

class HBsep(object):

    """
    i: data
    j: template
    k: class
    z: redshift
    f: filter
    """


    def __init__(self, data = None,
                 magcols = None,
                 errcols = None,
                 filters = None,
                 Nclasses = None,
                 limmaginfo = (None,None,None),
                 mismaginfo = None):

        self.data = self.read_data(self.data)
        self.process_data()

        self.Ndata    = self.data.shape[0]
        self.Nclasses = Nclasses

        self.filters = self.read_data(filters)
        

        


    def read_filters(self):
        

    def read_data(self,data):
        """
        Read in data file or array, do some moderate
        error checking.
        """
        if data==None:
            print '\nSpecification for data not given,'
            print 'please give a file location or 2D numpy array'
            assert False
        elif isinstance(data,np.ndarray):
            self.process_data()
        elif isinstance(data,str):
            try:
                # hdu number might be an issue
                f = pf.open(data)
                data = f[1].data
                f.close()
                self.shape_check()
            except:
                pass
            try:
                data = np.loadtxt(data)
            except:
                print '\nData file not read.'
                assert False

        self.shape_check(data,message)
        return data

    def shape_check(self,data,message):
        """
        Simple shape checks
        """
        assert len(data.shape)==2, \
            '\nData must be a 2D numpy array'
        assert data.shape[0]>data.shape[1], \ 
            '\nMust have more'
        assert np.mod(data.shape[1],2)==0, '\nNcolumn is not even,' + \
            'numpy array must have columns of mags then errs.'


    def process_data(self)
        """
        Turn data into mags, mag err, flux, flux err
        """
        self.Nfilter = self.data.shape[1] / 2
        if keepmags:
            self.mags = self.data[:,:Nfilter]
            self.merr = self.data[:,Nfilter:]

        for f in range(self.Nfilter):

            # Treat missing data
            # RF - need to vet
            if mismaginfo[0]!=None:
                ind = self.mags[:,f]==mismaginfo[f]
                print 'Found %d missing data in filter %d' % \
                    (len(self.mags[ind]),f)
                self.flux[ind,f] = 10.**(-0.4 * 20.0)  
                self.ferr[ind,f] = 0.4 * np.log(10.) * \
                    self.flux[ind,f] * 100. 
            else:
                print 'Assuming no missing data'

            #Treat undetected
            if limmaginfo[0]!=None:
                ind = self.mags[:,f]==limmaginfo[0,f]
                print 'Found %d undetected sources in filter %d' % \
                    (len(self.mags[ind],f))
                self.flux[ind,f] = 0.0
                self.ferr[ind,f] = 0.4 * np.log(10.) * \
                    10.**(-0.4 * limmaginfo[1,f]) * self.merr[ind,f]
            else:
                print 'Assuming no undetected sources'

            ind = (self.mags[:,f]>0.0) & (self.mags[:,f]<35.) & \
                (self.mags[:,f]!=mismaginfo[f]) & \
                (self.mags[:,f]!=limmaginfo[0,f])
            self.flux[ind,f] = 10.**(-0.4 * self.mags[ind,f])
                self.ferr[ind,f] = 0.4 * np.log(10.) * \
                    self.flux[ind,f] * self.merr[ind,f] 



            
