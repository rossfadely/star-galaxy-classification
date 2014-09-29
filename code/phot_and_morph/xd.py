"""
XD version originally by Jake Vanderplas/astroML:
https://github.com/astroML/astroML/blob/master/astroML/density_estimation/xdeconv.py

Modifications by Ross Fadely:
- 2014-09-23: port and erase astroML dependencies.
- 2014-09-25: convert class to functions for easy multiprocessing,
              only Estep seems to be useful across mutliple threads.
- 2014-09-29: add simple covarince regularization to avoid singular matricies.

Extreme deconvolution solver

This follows Bovy et al.
http://arxiv.org/pdf/0905.2979v2.pdf

Arbitrary mixing matrices R are not yet implemented: currently, this only
works with R = I.
"""

import multiprocessing
import numpy as np

from time import time
from sklearn.mixture import GMM
from utils import logsumexp, log_multivariate_gaussian

def XDGMM(X, Xcov, n_components, n_iter=100, tol=1E-5, Nthreads=1, R=None, 
          init_n_iter=10, w=None, verbose=False):
    """
    Extreme Deconvolution

    Fit an extreme deconvolution (XD) model to the data

    Parameters
    ----------
    n_components: integer
        number of gaussian components to fit to the data
    n_iter: integer (optional)
        number of EM iterations to perform (default=100)
    tol: float (optional)
        stopping criterion for EM iterations (default=1E-5)
    X:    array_like
          Input data. shape = (n_samples, n_features)
    Xcov: array_like
          Covariance of input data.  shape = (n_samples, n_features, n_features)
    R:    array_like
          (TODO: not implemented)
          Transformation matrix from underlying to observed data.  If
          unspecified, then it is assumed to be the identity matrix.
    w:    float or array_like
          if float - w * np.eye is added to V
          if vector - np.diag(w) is added to V
          if array - w is added to V 
    Notes
    -----
    This implementation follows Bovy et al. arXiv 0905.2979
    """
    model = xd_model(X.shape, n_components, n_iter, tol, w, Nthreads, verbose)

    if R is not None:
        raise NotImplementedError("mixing matrix R is not yet implemented")

    X = np.asarray(X)
    Xcov = np.asarray(Xcov)

    # assume full covariances of data
    assert Xcov.shape == (model.n_samples, model.n_features, model.n_features)

    # initialize components via a few steps of GMM
    # this doesn't take into account errors, but is a fast first-guess
    if model.V is None:
        t0 = time()
        gmm = GMM(model.n_components, n_iter=init_n_iter, covariance_type='full').fit(X)
        model.mu = gmm.means_
        model.alpha = gmm.weights_
        model.V = gmm.covars_

        if model.verbose:
            print 'Initalization done in %.2g sec' % (time() - t0)

    logL = model.logL(X, Xcov)

    for i in range(model.n_iter):
        t0 = time()
        model = _EMstep(model, X, Xcov)
        logL_next = model.logL(X, Xcov)
        t1 = time()

        if model.verbose:
            print "%i: log(L) = %.5g" % (i + 1, logL_next)
            print "    (%.2g sec)" % (t1 - t0)

        if logL_next < logL + model.tol:
            break
        logL = logL_next

    return model

def _EMstep(model, X, Xcov):
    """
    Perform the E-step (eq 16 of Bovy et al)
    """
    X = X[:, np.newaxis, :]
    Xcov = Xcov[:, np.newaxis, :, :]

    w_m = X - model.mu
    T = Xcov + model.V

    if model.Nthreads > 1:
        q, b, B = _Estep_multi(model, T, w_m, X)
    elif model.Nthreads == 1:
        q, b, B = _Estep((T, w_m, X, model, model.n_samples))
    else:
        assert False, 'Number of threads must be greater than 1.'

    return _Mstep(model, q, b, B)

def _Estep((T, w_m, X, model, n_samples)):
    """
    Compute the Estep for the given data chunk and current model
    """
    # compute inverse of each covariance matrix T
    Tshape = T.shape
    T = T.reshape([n_samples * model.n_components,
                   model.n_features, model.n_features])
    Tinv = np.array([np.linalg.inv(T[i])
                     for i in range(T.shape[0])]).reshape(Tshape)
    T = T.reshape(Tshape)

    # evaluate each mixture at each point
    N = np.exp(log_multivariate_gaussian(X, model.mu, T, Vinv=Tinv))

    # q
    q = (N * model.alpha) / np.dot(N, model.alpha)[:, None]
    
    # b
    tmp = np.sum(Tinv * w_m[:, :, np.newaxis, :], -1)
    b = model.mu + np.sum(model.V * tmp[:, :, np.newaxis, :], -1)

    # B
    tmp = np.sum(Tinv[:, :, :, :, np.newaxis]
                 * model.V[:, np.newaxis, :, :], -2)
    B = model.V - np.sum(model.V[:, :, :, np.newaxis]
                         * tmp[:, :, np.newaxis, :, :], -2)

    return (q, b, B)

def _Estep_multi(model, T, w_m, X):
    """
    Use multiple processes to compute Estep.
    """
    pool = multiprocessing.Pool(model.Nthreads)
    mapfn = pool.map
    Nchunk = np.ceil(1. / model.Nthreads * model.n_samples).astype(np.int)

    arglist = [None] * model.Nthreads
    for i in range(model.Nthreads):
        s = i * Nchunk
        e = s + Nchunk
        arglist[i] = (T[s:e], w_m[s:e], X[s:e], model, X[s:e].shape[0])

    results = list(mapfn(_Estep, [args for args in arglist]))

    q, b, B = results[0]
    for i in range(1, model.Nthreads):
        q = np.vstack((q, results[i][0]))
        b = np.vstack((b, results[i][1]))
        B = np.vstack((B, results[i][2]))

    pool.close()
    pool.terminate()
    pool.join()
    return q, b, B

def _Mstep(model, q, b, B):
    """
    M-step: compute alpha, mu, V, update to model
    """
    qj = q.sum(0)

    # update alpha
    model.alpha = qj / model.n_samples

    # update mu
    model.mu = np.sum(q[:, :, np.newaxis] * b, 0) / qj[:, np.newaxis]

    # update V
    m_b = model.mu - b
    tmp = m_b[:, :, np.newaxis, :] * m_b[:, :, :, np.newaxis]
    tmp += B
    tmp *= q[:, :, np.newaxis, np.newaxis]
    model.V = (tmp.sum(0) + model.w[np.newaxis, :, :]) / (qj[:, np.newaxis, np.newaxis] + 1)

    return model

class xd_model(object):
    """
    Class to store all things pertinent to the XD model. 
    """
    def __init__(self, xshape, n_components, n_iter, tol, w, Nthreads, verbose):
        self.n_samples = xshape[0]
        self.n_features = xshape[1]
        self.n_components = n_components
        self.n_iter = n_iter
        self.tol = tol
        self.Nthreads = Nthreads
        self.verbose = verbose

        self.V = None
        self.mu = None
        self.alpha = None

        # construct simple cov regularization term.  
        # regularize only along the diagonal, and same for each component.
        if type(w) == float:
            w = np.eye(self.n_features) * w
        elif type(w) == np.ndarray:
            if w.size == self.n_features:
                w = np.diag(w)
        else:
            w = np.diag(np.zeros(self.n_features))
        self.w = w

    def logprob_a(self, X, Xcov):
        """
        Evaluate the probability for a set of points

        Parameters
        ----------
        X: array_like
            Input data. shape = (n_samples, n_features)
        Xcov: array_like
            Covariance of input data.  shape = (n_samples, n_features, n_features)

        Returns
        -------
        p: ndarray
            Probabilities.  shape = (n_samples,)
        """
        X = np.asarray(X)
        Xcov = np.asarray(Xcov)
        n_samples, n_features = X.shape

        # assume full covariances of data
        assert Xcov.shape == (n_samples, n_features, n_features)

        X = X[:, np.newaxis, :]
        Xcov = Xcov[:, np.newaxis, :, :]
        T = Xcov + self.V

        return log_multivariate_gaussian(X, self.mu, T)

    def logL(self, X, Xcov):
        """
        Compute the log-likelihood of data given the model

        Parameters
        ----------
        X: array_like
            data, shape = (n_samples, n_features)
        Xcov: array_like
            covariances, shape = (n_samples, n_features, n_features)

        Returns
        -------
        logL : float
            log-likelihood
        """
        return np.sum(logsumexp(self.logprob_a(X, Xcov), -1))

    def sample(self, alpha, mu, V, size=1):
        shape = tuple(np.atleast_1d(size)) + (mu.shape[1],)
        npts = np.prod(size)

        alpha_cs = np.cumsum(alpha)
        r = np.atleast_1d(np.random.random(size))
        r.sort()

        ind = r.searchsorted(alpha_cs)
        ind = np.concatenate(([0], ind))
        if ind[-1] != size:
            ind[-1] = size

        draw = np.vstack([np.random.multivariate_normal(mu[i],
                                                        V[i],
                                                        (ind[i + 1] - ind[i],))
                          for i in range(len(alpha))])

        return draw.reshape(shape)

    def posterior(self, X, Xcov, sample=True):
        """
        Return the posterior mean and covariance given the xd model.
        """
        draw = None
        post_mu = np.zeros_like(X.shape[0], self.alpha.shape, X.shape[1])
        post_V = np.zeros_like(Xcov.shape[0], self.alpha.shape,
                               Xcov.shape[1], Xcov.shape[2])

        iV = np.zeros_like(self.V)
        for i in range(self.n_components):
            iV[i] = np.linalg.inv(self.V[i])

        for i in range(Xcov.shape[0]):
            Xicov = np.linalg.inv(Xcov[i])
            for j in range(self.n_components):
                post_V[i, j] = np.linalg.inv(iV[j] + Xicov)
                post_mu[i, j] = np.dot(iV[j], self.mu[j])
                post_mu[i, j] += np.dot(Xicov, X[i])
                post_mu[i, j] = np.dot(post_V[i, j], post_mu[i, j])

if __name__ == '__main__':

    from matplotlib import use; use('Agg')
    from utils import fetch_mixed_epoch
    
    import os
    import matplotlib.pyplot as pl

    epoch = 10
    single, coadd = fetch_mixed_epoch(epoch)

    psfmags = np.vstack([single['psfmag_' + f] - 
                         single['extinction_' + f] for f in 'ugriz']).T
    psfmagerrs = np.vstack([single['psfmagerr_' + f] for f in 'ugriz']).T
    modelmags = np.vstack([single['modelmag_' + f] -
                           single['extinction_' + f]for f in 'ugriz']).T
    modelmagerrs = np.vstack([single['modelmagerr_' + f] for f in 'ugriz']).T
    
    X = np.hstack((psfmags, modelmags))
    Xerr = np.hstack((psfmagerrs, modelmagerrs))

    W = np.zeros((10, X.shape[1]))
    W[0, 2] = 1.
    for i in range(1, 6):
        W[i, i - 1] = 1
        W[i, i + 4] = -1
    for i in range(6, 10):
        W[i, i - 1] = 1
        W[i, i] = -1

    X = np.dot(X, W.T)
    Xcov = np.zeros(Xerr.shape + Xerr.shape[-1:])
    Xcov[:, range(Xerr.shape[1]), range(Xerr.shape[1])] = Xerr ** 2
    Xcov = np.tensordot(np.dot(Xcov, W.T), W, (-2, -1))

    if False:
        w = np.ones(X.shape[1]) * np.inf
        for i in range(X.shape[0]):
            w = np.minimum(w, np.diag(Xcov[i]))
        w /= 10.
    else:
        w = None

    xd = XDGMM(X, Xcov, 32, n_iter=40, w=w, Nthreads=10, verbose=True)

    import cPickle
    f=open(os.environ['sgdata'] + 'xd_32_%d.pkl' % epoch, 'wb')
    cPickle.dump(xd, f)
    f.close()
