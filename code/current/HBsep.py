
import numpy as np
import pyfits as pf
import matplotlib.pyplot as pl

class HBsep(object):

    """
    i: data
    j: template
    k: class
    z: redshift
    f: filter
    m: misc
    """


    def __init__(self, data,
                 magcols = None,
                 errcols = None,
                 keepmags = False,
                 filters = None,
                 Nclasses = None,
                 limmaginfo = (None,None,None),
                 mismaginfo = None):

        self.mismaginfo = mismaginfo
        self.limmaginfo = limmaginfo

        print 'here'
        data = self.read_data(data,'Data Error:')
        self.Ndata   = data.shape[0]
        self.Nfilter = data.shape[1] / 2
        self.mags = data[:,:self.Nfilter]
        self.merr = data[:,self.Nfilter:]
        self.flux = self.mags * 0.0
        self.ferr = self.merr * 0.0

        self.process_data()
        #self.Nclasses = Nclasses

        #self.filters = self.read_data(filters)

    def read_data(self,data,item):
        """
        Read in data file or array, do some moderate
        error checking.
        """
        if isinstance(data,np.ndarray):
            pass
        elif isinstance(data,str):
            try:
                # hdu number might be an issue
                f = pf.open(data)
                tdata = f[1].data
                f.close()
                data = np.zeros((len(tdata.field(0)),
                                 len(tdata[0])))
                for m in range(data.shape[1]):
                    data[:,m] = tdata.field(m)                
            except:
                try:
                    data = np.loadtxt(data)
                except:
                    print '\nData file not read.'
                    assert False

        self.shape_check(data,item)
        return data

    def shape_check(self,data,item):
        """
        Simple shape checks
        """
        assert len(data.shape)==2, '\n\n'+item+' Must be a 2D numpy array'
        assert data.shape[0]>data.shape[1], '\n\n'+item+' Must have more rows than columns'
        assert np.mod(data.shape[1],2)==0, '\n\n'+item+' Ncolumn is not even, what gives?'


    def process_data(self):
        """
        Turn data into mags, mag err, flux, flux err
        """
        for f in range(self.Nfilter):

            # Treat missing data
            # RF - need to vet
            #if self.mismaginfo[0]!=None:
            #    ind = self.mags[:,f]==self.mismaginfo[f]
            #    print 'Found %d missing data in filter %d' % \
            #        (len(self.mags[ind]),f)
            #    self.flux[ind,f] = 10.**(-0.4 * 20.0)  
            #    self.ferr[ind,f] = 0.4 * np.log(10.) * \
            #        self.flux[ind,f] * 100. 
            #else:
            #    print 'Assuming no missing data in filter %d' % f

            #Treat undetected
            #if self.limmaginfo[0]!=None:
            #    ind = self.mags[:,f]==limmaginfo[0,f]
            #    print 'Found %d undetected sources in filter %d' % \
            #        (len(self.mags[ind],f))
            #    self.flux[ind,f] = 0.0
            #    self.ferr[ind,f] = 0.4 * np.log(10.) * 10.**(-0.4 * self.limmaginfo[1,f]) * self.merr[ind,f]
            #else:
            #    print 'Assuming no undetected sources in filter %d' % f

            #Process rest
            #ind = (self.mags[:,f]>0.0) & (self.mags[:,f]<35.) & \
            #    (self.mags[:,f]!=self.mismaginfo[f])# & \
            #    #(self.mags[:,f]!=self.limmaginfo[0,f])
            #self.flux[ind,f] = 10.**(-0.4 * self.mags[ind,f])
            #self.ferr[ind,f] = 0.4 * np.log(10.) * self.flux[ind,f] * self.merr[ind,f] 
            pass
        self.flux = 10.**(-0.4 * self.mags)
        self.ferr = 0.4 * np.log(10.)
        print self.flux
            
if __name__=='__main__':
    HBsep('cat_w_pt_all.fits')
