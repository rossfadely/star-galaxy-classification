# 
# Various routines, mostly for handing data.
#
# Author - Ross Fadely
#
import os
import numpy as np
import pyfits as pf

def fetch_epoch(epoch, kind, verbose=False):
    """
    Return the single epoch and the matched coadded data.
    """
    assert kind in ['stars', 'gals']
    ddir = os.environ['sgdata']

    # single epoch
    f = pf.open(ddir + 's82single_%s_%d.fits' % (kind, epoch))
    s = f[1].data
    f.close()
    N = s.field(0).size

    try:
        f = pf.open(ddir + 's82coadd_%s_%d.fits' % (kind, epoch))
        c = f[1].data
        f.close()

    except:
        # master
        f = pf.open(ddir + 's82coadd30k_%s_rfadely.fit' % kind)
        c = f[1].data
        f.close()

        # find the corresponding coadds
        inds = np.zeros(N, dtype=np.int)
        ind = 0
        for i in range(N):
            if verbose:
                if i % 200 == 0:
                    print 'searching', i
            coadd_objid = s.field('coadd_objid')[i]
            search = True
            while search:
                if c.field('objid')[ind] == coadd_objid:
                    inds[i] = ind
                    search = False
                ind += 1

        c = c[inds]
        if False:
            # paranoid check
            for i in range(N):
                st = '%d %d' % (c[i].field('objid'), s[i].field('coadd_objid'))
                assert c[i].field('objid') == s[i].field('coadd_objid'), st

        dt = {'E':np.float32, 'K':np.int64, 'D':np.float64, 'I':np.int16, 'K':np.int64}
        cols = []
        for i in range(len(s[0]) - 1):
            n = s._coldefs.names[i]
            f = s._coldefs.formats[i]
            cols.append(pf.Column(name=n, format=f, array=c.field(n).astype(dt[f])))

        tbhdu = pf.new_table(pf.ColDefs(cols))
        tbhdu.writeto(ddir + 's82coadd_%s_%d.fits' % (kind, epoch), clobber=True)
        c = tbhdu.data
    
    return c

if __name__ == '__main__':
    import time
    t=time.time()
    fetch_epoch(12, 'stars', True)
    print time.time()-t
