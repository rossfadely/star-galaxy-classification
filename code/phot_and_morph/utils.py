# 
# Various routines, mostly for handing data.
#
# Author - Ross Fadely
#
import os
import multiprocessing
import numpy as np
import pyfits as pf

def log_multivariate_gaussian_Nthreads(x, mu, V, Nthreads=1):
    """
    Use multiprocessing to calculate log likelihoods.
    """
    n_samples = x.shape[0]

    pool = multiprocessing.Pool(Nthreads)
    mapfn = pool.map
    Nchunk = np.ceil(1. / Nthreads * n_samples).astype(np.int)

    if len(V.shape) > 3:
        perdatum = True
    else:
        perdatum = False

    arglist = [None] * Nthreads
    for i in range(Nthreads):
        s = i * Nchunk
        e = s + Nchunk
        if perdatum:
            arglist[i] = (x[s:e], mu, V[s:e])
        else:
            arglist[i] = (x[s:e], mu, V)

    result = list(mapfn(lmg, [args for args in arglist]))
    
    logls = result[0]
    for i in range(1, Nthreads):
       logls = np.vstack((logls, result[i]))
       
    pool.close()
    pool.terminate()
    pool.join()
    return logls

def lmg(args):
    return log_multivariate_gaussian(*args)

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
        Vchol = np.array([np.linalg.cholesky(V[i])
                          for i in range(V.shape[0])])

        # we may be more efficient by using scipy.np.linalg.solve_triangular
        # with each cholesky decomposition
        VcholI = np.array([np.linalg.inv(Vchol[i])
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
            Vinv = np.array([np.linalg.inv(V[i])
                             for i in range(V.shape[0])]).reshape(Vshape)
        else:
            assert Vinv.shape == Vshape

        logdet = np.log(np.array([np.linalg.det(V[i])
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
        print 'Matched fits for coadd doesn\'t exist, building...'

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

def fetch_prepped_s82data(epoch, fgal=0.5, features=['psf_mag', 'model_colors',
                                                     'psf_minus_model'],
                          filters=['r', 'ur gr ri rz', 'r'], use_single=True):
    """
    Construct data matrix and cov.
    """
    single, coadd = fetch_mixed_epoch(epoch, fgal)

    if use_single:
        d = single
    else:
        d = coadd
    return prep_data(d, features, filters)

def fetch_prepped_dr10data(N, fgal=0.5, features=['psf_mag', 'model_colors',
                                                  'psf_minus_model'],
                           filters=['r', 'ur gr ri rz', 'r'],
                           seed=1234):
    """
    Prepare SDSS DR10 data to run XD.
    """
    np.random.seed(seed)
    ddir = os.environ['sgdata']
    f = pf.open(ddir + 'dr10_30k_stars.fits')
    d = f[1].data
    f.close()
    Xstar, Xstarcov = prep_data(d, features, filters)
    f = pf.open(ddir + 'dr10_30k_gals.fits')
    d = f[1].data
    f.close()
    Xgal, Xgalcov = prep_data(d, features, filters)

    Nmax = np.minimum(Xstar.shape[0], Xgal.shape[0])
    assert N <= Nmax, 'Not enough data for request'

    Ngal = np.round(N * fgal).astype(np.int)
    Nstar = N - Ngal
    ind = np.random.permutation(Xgal.shape[0])[:Ngal]
    Xgal = Xgal[ind]
    Xgalcov = Xgalcov[ind]
    ind = np.random.permutation(Xstar.shape[0])[:Nstar]
    Xstar = Xstar[ind]
    Xstarcov = Xstarcov[ind]

    X = np.vstack((Xgal, Xstar))
    Xcov = np.vstack((Xgalcov, Xstarcov))
    ind = np.random.permutation(X.shape[0])
    X = X[ind]
    Xcov = Xcov[ind]
    return X, Xcov

def make_W_matrix(features, filters, odim):
    """
    Construct the mixing matrix for the set of features.
    """
    ref = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4}
    W = np.zeros((1, odim))
    idx = 0
    for i, feature in enumerate(features):
        # begin spaghetti code
        if 'psf' in feature:
            ind = 0
        else:
            ind = 5

        if 'mag' in feature:
            for f in filters[i]:
                W = np.vstack((W, np.zeros((1, odim))))
                W[idx, ind + ref[f]] = 1.
                idx += 1

        if 'colors' in feature:
            filts = filters[i].split()
            for f in filts:
                W = np.vstack((W, np.zeros((1, odim))))
                W[idx, ind + ref[f[0]]] = 1.
                W[idx, ind + ref[f[1]]] = -1.
                idx += 1
        
        if 'minus' in feature:
            for f in filters[i]:
                W = np.vstack((W, np.zeros((1, odim))))
                W[idx, ref[f]] = 1.
                W[idx, 5 + ref[f]] = -1.
                idx += 1
            
    return W[:-1]

def prep_data(d, features, filters=None):
    """
    Return the prepared data.
    """
    if filters is None:
        filters = ['ugriz' for i in range(len(features))]

    psfmags = np.vstack([d['psfmag_' + f] -
                         d['extinction_' + f] for f in 'ugriz']).T
    psfmagerrs = np.vstack([d['psfmagerr_' + f] for f in 'ugriz']).T
    modelmags = np.vstack([d['modelmag_' + f] -
                           d['extinction_' + f]for f in 'ugriz']).T
    modelmagerrs = np.vstack([d['modelmagerr_' + f] for f in 'ugriz']).T

    X = np.hstack((psfmags, modelmags))
    Xerr = np.hstack((psfmagerrs, modelmagerrs))

    W = make_W_matrix(features, filters, X.shape[1])

    X = np.dot(X, W.T)
    Xcov = np.zeros(Xerr.shape + Xerr.shape[-1:])
    Xcov[:, range(Xerr.shape[1]), range(Xerr.shape[1])] = Xerr ** 2
    Xcov = np.tensordot(np.dot(Xcov, W.T), W, (-2, -1))

    return X, Xcov

if __name__ == '__main__':
    import time
    t=time.time()
    fetch_epoch(12, 'stars', True)
    print time.time()-t
