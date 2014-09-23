# 
# Various routines, mostly for handing data.
#
# Author - Ross Fadely
#
import os
import numpy as np
import pyfits as pf

def log_multivariate_gaussian(x, mu, V, Vinv=None, method=1):
    """
    Swiped from astroML:
    https://github.com/astroML/astroML/blob/master/astroML/utils.py

    Evaluate a multivariate gaussian N(x|mu, V)

    This allows for multiple evaluations at once, using array broadcasting

    Parameters
    ----------
    x: array_like
        points, shape[-1] = n_features

    mu: array_like
        centers, shape[-1] = n_features

    V: array_like
        covariances, shape[-2:] = (n_features, n_features)

    Vinv: array_like or None
        pre-computed inverses of V: should have the same shape as V

    method: integer, optional
        method = 0: use cholesky decompositions of V
        method = 1: use explicit inverse of V

    Returns
    -------
    values: ndarray
        shape = broadcast(x.shape[:-1], mu.shape[:-1], V.shape[:-2])

    Examples
    --------

    >>> x = [1, 2]
    >>> mu = [0, 0]
    >>> V = [[2, 1], [1, 2]]
    >>> log_multivariate_gaussian(x, mu, V)
    -3.3871832107434003
    """
    x = np.asarray(x, dtype=float)
    mu = np.asarray(mu, dtype=float)
    V = np.asarray(V, dtype=float)

    ndim = x.shape[-1]
    x_mu = x - mu

    if V.shape[-2:] != (ndim, ndim):
        raise ValueError("Shape of (x-mu) and V do not match")

    Vshape = V.shape
    V = V.reshape([-1, ndim, ndim])

    if Vinv is not None:
        assert Vinv.shape == Vshape
        method = 1

    if method == 0:
        Vchol = np.array([linalg.cholesky(V[i], lower=True)
                          for i in range(V.shape[0])])

        # we may be more efficient by using scipy.linalg.solve_triangular
        # with each cholesky decomposition
        VcholI = np.array([linalg.inv(Vchol[i])
                          for i in range(V.shape[0])])
        logdet = np.array([2 * np.sum(np.log(np.diagonal(Vchol[i])))
                           for i in range(V.shape[0])])

        VcholI = VcholI.reshape(Vshape)
        logdet = logdet.reshape(Vshape[:-2])

        VcIx = np.sum(VcholI * x_mu.reshape(x_mu.shape[:-1]
                                            + (1,) + x_mu.shape[-1:]), -1)
        xVIx = np.sum(VcIx ** 2, -1)

    elif method == 1:
        if Vinv is None:
            Vinv = np.array([linalg.inv(V[i])
                             for i in range(V.shape[0])]).reshape(Vshape)
        else:
            assert Vinv.shape == Vshape

        logdet = np.log(np.array([linalg.det(V[i])
                                  for i in range(V.shape[0])]))
        logdet = logdet.reshape(Vshape[:-2])

        xVI = np.sum(x_mu.reshape(x_mu.shape + (1,)) * Vinv, -2)
        xVIx = np.sum(xVI * x_mu, -1)

    else:
        raise ValueError("unrecognized method %s" % method)

    return -0.5 * ndim * np.log(2 * np.pi) - 0.5 * (logdet + xVIx)

def logsumexp(arr, axis=None):
    """
    Swiped from astroML:
    https://github.com/astroML/astroML/blob/master/astroML/utils.py
    
    Computes the sum of arr assuming arr is in the log domain.

    Returns log(sum(exp(arr))) while minimizing the possibility of
    over/underflow.

    Examples
    --------

    >>> import numpy as np
    >>> a = np.arange(10)
    >>> np.log(np.sum(np.exp(a)))
    9.4586297444267107
    >>> logsumexp(a)
    9.4586297444267107
    """
    # if axis is specified, roll axis to 0 so that broadcasting works below
    if axis is not None:
        arr = np.rollaxis(arr, axis)
        axis = 0

    # Use the max to normalize, as with the log this is what accumulates
    # the fewest errors
    vmax = arr.max(axis=axis)
    out = np.log(np.sum(np.exp(arr - vmax), axis=axis))
    out += vmax

    return out

def fetch_mixed_epoch(epoch, gal_frac=0.5, shuffle=True, seed=12345):
    """
    Make a combined sample of stars and galaxies for a given epoch.  Return single epoch and 
    coadd data.
    """
    # fetch
    ss, sc = fetch_epoch(epoch, 'stars')
    gs, gc = fetch_epoch(epoch, 'gals')

    # size
    N = np.minimum(ss.field(0).size, gs.field(0).size)
    Ngal = np.round(gal_frac * N)
    Nstar = N - Ngal

    # shuffle:
    if shuffle:
        np.random.seed(seed)
        ind = np.random.permutation(N).astype(np.int)
    else:
        ind = np.arange(N, dtype=np.int)

    # cut
    ss, sc = ss[:Nstar], sc[:Nstar]
    gs, gc = gs[:Ngal], gc[:Ngal]

    # build
    os = {}
    oc = {}
    for k in ss._coldefs.names:
        os[k] = np.append(ss.field(k), gs.field(k))[ind]
        if k != 'coadd_objid':
            oc[k] = np.append(sc.field(k), gc.field(k))[ind]

    return os, oc

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
    
    return s, c

if __name__ == '__main__':
    import time
    t=time.time()
    fetch_epoch(12, 'stars', True)
    print time.time()-t
