#
# Plotting routines, mostly for exploration.  Paper figures are generated with `fig_x`
#
# Author - Ross Fadely
#
from matplotlib import use; use('Agg')
from utils import fetch_mixed_epoch, fetch_prepped_s82data
from xd import xd_model

import os
import cPickle
import numpy as np
import matplotlib.pyplot as pl

def cmd_plot(cmd1, cmd2, titles, fname):
    """
    Plot the cmd for two different datasets.
    """
    xmn, xmx = cmd1[:, 0].min(), cmd1[:, 0].max()
    ymn, ymx = cmd1[:, 1].max(), cmd1[:, 1].min()
    a, fs = 0.75, 3
    f = pl.figure(figsize=(2 * fs, fs))
    pl.subplot(121)
    pl.scatter(cmd1[:, 0], cmd1[:, 1], marker='o', color='k', alpha=a)
    pl.xlim(xmn, xmx)
    pl.ylim(ymn, ymx)
    pl.title(titles[0])
    pl.subplot(122)
    pl.scatter(cmd2[:, 0], cmd2[:, 1], marker='o', color='k', alpha=a)
    pl.xlim(xmn, xmx)
    pl.ylim(ymn, ymx)
    pl.title(titles[1])
    f.savefig(fname)

def s82_compare(data, post, s82, name, lim=np.array([0.025, 0.975]),
                plotdir='../../plots/phot_and_morph/'):
    """
    Plot the ugriz comparison to Stripe 82 data.
    """
    bands = 'ugriz'
    fs = 5
    mn, mx = 18., 24
    f = pl.figure(figsize=(5 * fs, 2 * fs))
    for i in range(5):
        pl.subplot(2, 5, i + 1)
        vs = np.sort(data[:, i])
        #mn, mx = vs[lim[0] * vs.size], vs[lim[1] * vs.size] 
        pl.plot(data[:, i], s82[:, i], 'k.', alpha=0.5)
        pl.plot([mn, mx], [mn, mx], 'r', lw=2)
        pl.xlim(mn, mx)
        pl.ylim(mn, mx)
        pl.title(bands[i], fontsize=20)
        pl.xlabel('Single Epoch')
        if i == 0:
            pl.ylabel('Coadd')
    for i in range(5):
        pl.subplot(2, 5, i + 6)
        vs = np.sort(data[:, i])
        #mn, mx = vs[lim[0] * vs.size], vs[lim[1] * vs.size]
        pl.plot(post[:, i], s82[:, i], 'k.', alpha=0.5)
        pl.plot([mn, mx], [mn, mx], 'r', lw=2)
        pl.xlim(mn, mx)
        pl.ylim(mn, mx)
        pl.title(bands[i], fontsize=20)
        pl.xlabel('XD Posterior')
        if i == 0:
            pl.ylabel('Coadd')
    f.savefig(plotdir + name)

def plot_misclass(epoch=10, plotdir='../../plots/phot_and_morph/'):
    """
    Just plot the single epoch `galaxies` that turned out to be 
    stars in the coadd - its weird that this happens at bright mags.
    """
    single, coadd = fetch_mixed_epoch(epoch)

    mags = coadd['psfmag_r']
    singlepmm = single['psfmag_r'] - single['modelmag_r']
    coaddpmm = coadd['psfmag_r'] - coadd['modelmag_r']

    ind = (single['type'] ==  3) & (coaddpmm < 0.03) & (coaddpmm > -0.03)

    f = pl.figure()
    pl.plot(mags[ind], singlepmm[ind], '.k', alpha=0.3)
    f.savefig(plotdir + 'misclassed_gals.png')

def one_epoch_class_check(epoch=10, plotdir='../../plots/phot_and_morph/'):
    """
    Check out the change in classification from the single epoch to the coadd.
    """
    single, coadd = fetch_mixed_epoch(epoch)

    ind = np.random.permutation(single['objid'].size)[0]
    assert single['coadd_objid'][ind] == coadd['objid'][ind], 'Object mismatch'

    mags = coadd['psfmag_r']
    singlepmm = single['psfmag_r'] - single['modelmag_r']
    coaddpmm = coadd['psfmag_r'] - coadd['modelmag_r']

    fs = 5
    alpha = 0.3
    colors = ['#CC6600', '#990099']
    titles = ['single epoch ', 'coadd ']
    classes = [single['type'], coadd['type']]
    xlabel = 'Coadd psfmag r'
    ylabel = 'psf - model r'
    xlim = (18, 22)
    ylim = (-0.1, 0.5)
    
    f = pl.figure(figsize=(4 * fs, 2 * fs))

    pl.subplots_adjust(wspace=0.25, right=0.8)

    xticks = []
    for i in range(10):
        if i % 2 == 0:
            xticks.append('%0.0f' % (18 + 0.5 * i))
        else:
            xticks.append('')

    types = [3, 6]
    crit = [0.145, 0.03]
    for i in range(2):
        for j in range(2):
            ind = classes[i] == types[j]

            ax = pl.subplot(2, 4, 4 * i + j * 2 + 1)
            pl.plot(mags[ind], singlepmm[ind], '.', alpha=alpha, color=colors[j])
            pl.plot([18, 22], [crit[0], crit[0]], 'k', lw=2.5)
            pl.plot([18, 22], [-crit[0], -crit[0]], 'k', lw=2.5)
            pl.xlabel(xlabel)
            pl.ylabel(titles[0] + ylabel)
            pl.title(titles[i] + 'type')
            pl.xlim(xlim)
            pl.ylim(ylim)
            ax.set_xticklabels(xticks)

            ax = pl.subplot(2, 4, 4 * i + j * 2 + 2)
            pl.plot(mags[ind], coaddpmm[ind], '.', alpha=alpha, color=colors[j])
            pl.plot([18, 22], [crit[1], crit[1]], 'k', lw=2.5)
            pl.plot([18, 22], [-crit[1], -crit[1]], 'k', lw=2.5)
            pl.xlabel(xlabel)
            pl.ylabel(titles[1] + ylabel)
            pl.title(titles[i] + 'type')
            pl.xlim(xlim)
            pl.ylim(ylim)
            ax.set_xticklabels(xticks)

    rfnotes = 'RF:\n - Why are there single epoch "stars"\n so far above 0.145? '
    rfnotes += 'This seems to be\n true regardless of "single epoch type"\n or '
    rfnotes += '"coadd type".  Is the single epoch\n data crummier than normal SDSS?\n\n'
    rfnotes += '- It is plain to see that "type" from single\n epoch data is calling some '
    rfnotes += 'small\n fraction of stars "galaxies" and vice versa.\n\n'
    rfnotes += '- The coadd classification clears things up,\n but there must still be a '
    rfnotes += 'fair bit of label\n error.  "Galaxies" are piling up near 0.03\n in the coadd '
    rfnotes += 'and the dist. of stars is\n asymmetric at the faint end.\n\n'
    rfnotes += '- Note the downward trend of the stellar\n locus in the single epoch data.  '
    rfnotes += 'This bias\n is likely arises from the noise model used\n in the psf fitting.  '
    rfnotes += 'This is also seen in the\n coadd if you go to faint depths.  '
    rfnotes += 'This\n motivates not forcing axis alignment in all\n XD "star" clusters.'
    pl.text(1.1, 2.1, rfnotes, transform=ax.transAxes, va='top')

    f.savefig(plotdir + 'pminusm_changes_%d.png' % epoch)

def psf_minus_model_hists(epoch, model, features, filters, fgal=0.5,
                          idx=-1):
    """
    Plot the histogram of psf-model for the data and denoised data.
    """
    fs = 5
    w = 0.1
    nb = 50
    mags = [19.5, 20.5, 21.5]

    Xsingle, Xsinglecov = fetch_prepped_s82data(epoch, fgal, features, filters)
    Xcoadd, Xcoaddcov = fetch_prepped_s82data(epoch, fgal, features, filters,
                                              use_single=False)

    f = pl.figure(figsize=(3 * fs, fs))
    for i in range(len(mags)):
        ind = (Xsingle[:, 0] > mags[i] - w) & (Xsingle[:, 0] < mags[i] + w)
        #ind = ind & (Xsingle[:, idx] < 0.3)
        ind = ind & (Xcoadd[:, idx] < 0.03)
        a, m, v = model.posterior(Xsingle[ind], Xsinglecov[ind])
        posts = np.zeros_like(Xsingle[ind])
        for j in range(Xsingle[ind].shape[0]):
            posts[j] = np.median(model.sample(a[j], m[j], v[j], size=1000),
                                 axis=0)

        print i, mags[i]
        pl.subplot(1, 3, i + 1)
        pl.hist(Xsingle[ind, idx], nb, alpha=0.3, color='k')
        pl.hist(posts[:, idx], nb, alpha=0.3, color='r')
        pl.title('$r=%0.1f$' % mags[i])
        pl.xlim(-0.2, 0.2)
        pl.xlabel('psf - model')
    f.savefig('../../plots/phot_and_morph/foo.png')

def plot_contours_and_data(epoch, model, features, filters, fgal=0.5, idx=-1):
    """
    Plot the data and the contours for stars and galaxies.
    """
    from astroML.plotting.tools import draw_ellipse

    Xsingle, Xsinglecov = fetch_prepped_s82data(epoch, fgal, features, filters)
    Xcoadd, Xcoaddcov = fetch_prepped_s82data(epoch, fgal, features, filters,
                                              use_single=False)

    sind = Xcoadd[:, idx] < 0.03
    gind = Xcoadd[:, idx] > 0.03

    fs = 5
    ms = 1
    f = pl.figure(figsize=(3 * fs, 2 * fs))
    Nstar = len(np.where(model.fixed_means[:, idx] != np.inf)[0])
    pl.subplot(231)
    idx = [0, -1]
    for i in range(Nstar):
        print i, model.V[i, idx][:, idx]
        draw_ellipse(model.mu[i, idx], model.V[i, idx][:, idx], scales=[2],
                     ec='k', fc='gray', alpha=0.2)
    pl.plot(Xsingle[sind][::10, idx[0]], Xsingle[sind][::10, idx[1]], '.k',
            ms=ms)
    pl.xlim(18,22)
    pl.ylim(-0.1, 0.5)
    pl.subplot(232)
    idx = [2, 1]
    for i in range(Nstar):
        print i, model.V[i, idx][:, idx]
        draw_ellipse(model.mu[i, idx], model.V[i, idx][:, idx], scales=[2],
                     ec='k', fc='gray', alpha=0.2)
    pl.plot(Xsingle[sind][::10, idx[0]], Xsingle[sind][::10, idx[1]], '.k',
            ms=ms)
    pl.xlim(-2, 3)
    pl.ylim(-1, 6)
    pl.subplot(233)
    idx = [3, 4]
    for i in range(Nstar):
        draw_ellipse(model.mu[i, idx], model.V[i, idx][:, idx], scales=[2],
                     ec='k', fc='gray', alpha=0.2)
    pl.plot(Xsingle[sind][::10, idx[0]], Xsingle[sind][::10, idx[1]], '.k',
            ms=ms)
    pl.xlim(-2, 3)
    pl.ylim(-1, 3)
    pl.subplot(234)
    idx = [0, -1]
    for i in range(Nstar, model.n_components):
        print i, model.V[i, idx][:, idx]
        draw_ellipse(model.mu[i, idx], model.V[i, idx][:, idx], scales=[2],
                     ec='k', fc='gray', alpha=0.2)
    pl.plot(Xsingle[gind][::10, idx[0]], Xsingle[gind][::10, idx[1]], '.k',
            ms=ms)
    pl.xlim(18,22)
    pl.ylim(-0.1, 0.5)
    pl.subplot(235)
    idx = [2, 1]
    for i in range(Nstar, model.n_components):
        print i, model.V[i, idx][:, idx]
        draw_ellipse(model.mu[i, idx], model.V[i, idx][:, idx], scales=[2],
                     ec='k', fc='gray', alpha=0.2)
    pl.plot(Xsingle[gind][::10, idx[0]], Xsingle[gind][::10, idx[1]], '.k',
            ms=ms)
    pl.xlim(-2, 3)
    pl.ylim(-1, 6)
    pl.subplot(236)
    idx = [3, 4]
    for i in range(Nstar, model.n_components):
        draw_ellipse(model.mu[i, idx], model.V[i, idx][:, idx], scales=[2],
                     ec='k', fc='gray', alpha=0.2)
    pl.plot(Xsingle[gind][::10, idx[0]], Xsingle[gind][::10, idx[1]], '.k',
            ms=ms)
    pl.xlim(-2, 3)
    pl.ylim(-1, 3)
    f.savefig('../../plots/phot_and_morph/foo.png')

if __name__ == '__main__':
    #one_epoch_class_check()
    #plot_misclass()

    epoch = 3
    N = 30000
    Nr = epoch
    K = 32
    n_iter = 128
    Nstar = 20
    fixed_inds = [-1]
    data = 's82'
    factor = 1000.
    features = ['psf_mag', 'model_colors', 'psf_minus_model']
    filters = ['r', 'ur gr ri rz', 'r']
    message = 'pm_mc_pmm_r_all_r'
    fname = 'xdmodel_%s_%d_%d_%d_%d_%s.pkl' % (data, Nr, K, n_iter, Nstar,
                                               message)

    f = open(os.environ['sgdata'] + fname, 'rb')
    model = cPickle.load(f)
    f.close()
 
    if True:
        plot_contours_and_data(epoch, model, features, filters, idx=-1)

    if False:
        psf_minus_model_hists(epoch, model, features, filters, idx=-1)

    if False:
        epoch = 10
        Xsingle, Xsinglecov = fetch_prepped_s82data(epoch)
        Xcoadd, Xcoaddcov = fetch_prepped_s82data(epoch, use_single=False)

        N = 5000
        Xsingle = Xsingle[:N]
        Xsinglecov = Xsinglecov[:N]
        Xcoadd = Xcoadd[:N]
        Xcoaddcov = Xcoaddcov[:N]
        
        f = open(os.environ['sgdata'] + 'xd_w_32_10_2.pkl')
        model = cPickle.load(f)
        f.close()
        
        post_alpha, post_mu, post_V = model.posterior(Xsingle, Xsinglecov)
        posts = np.zeros_like(Xsingle)
        for i in range(N):
            posts[i] = np.median(model.sample(post_alpha[i], post_mu[i], post_V[i],
                                              size=1000), axis=0)

        Xsingle[:, 0] -= Xsingle[:, 3]  
        single = np.vstack((Xsingle[:, 6] + Xsingle[:, 0], 
                            Xsingle[:, 7] + Xsingle[:, 0],
                            Xsingle[:, 0],
                            Xsingle[:, 8] + Xsingle[:, 0],
                            Xsingle[:, 9] + Xsingle[:, 0])).T

        Xcoadd[:, 0] -= Xcoadd[:, 3]  
        coadd = np.vstack((Xcoadd[:, 6] + Xcoadd[:, 0], 
                            Xcoadd[:, 7] + Xcoadd[:, 0],
                            Xcoadd[:, 0],
                            Xcoadd[:, 8] + Xcoadd[:, 0],
                            Xcoadd[:, 9] + Xcoadd[:, 0])).T

        posts[:, 0] -= posts[:, 3]
        post = np.vstack((posts[:, 6] + posts[:, 0], 
                          posts[:, 7] + posts[:, 0],
                          posts[:, 0],
                          posts[:, 8] + posts[:, 0],
                          posts[:, 9] + posts[:, 0])).T

        s82_compare(single, post, coadd, 'foo.png', lim=np.array([0.025, 0.975]),
                    plotdir='../../plots/phot_and_morph/')
        """
        f=pl.figure(figsize=(10,5))
        pl.subplot(121)
        pl.plot(Xsingle[:, 0], Xsingle[:, 2] - Xcoadd[:, 2], '.k')
        pl.ylim(-.2,.2)
        pl.subplot(122)
        pl.plot(Xsingle[:, 0], posts[:, 2] - Xcoadd[:, 2], '.k')
        pl.ylim(-.2,.2)
        f.savefig('../../plots/phot_and_morph/foo.png')
        """

    if False:
        try:
            f = pf.__version__
        except:
            import pyfits as pf
        f = pf.open(os.environ['sgdata'] + 'seg3_1.2.fits')
        d = f[1].data
        f.close()
        ind = d.field('type') == 6
        d = d[ind]

        cmd = np.zeros((len(d.field(0)), 2))
        g = d.field('psfmag_g') - d.field('extinction_g')
        r = d.field('psfmag_r') - d.field('extinction_r')
        cmd[:, 0] = g - r
        cmd[:, 1] = r

        f = open(os.environ['sgdata'] + 'xdmodel_dr10_30000_32_128.pkl')
        model = cPickle.load(f)
        f.close()

        from utils import prep_data
        X, Xcov = prep_data(d, ['psf_minus_model', 'colors'])
        post_alpha, post_mu, post_V = model.posterior(X, Xcov)
        posts = np.zeros_like(X)
        for i in range(X.shape[0]):
            posts[i] = np.median(model.sample(post_alpha[i], post_mu[i], post_V[i],
                                              size=1000), axis=0)
            if post_alpha[i].sum() != post_alpha[i].sum():
                posts[i] = np.zeros_like(X[0])

        r = posts[:, 0]
        gmr = posts[:, 7]
        pcmd = np.zeros((X.shape[0], 2))
        pcmd[:, 0] = posts[:, 7]
        pcmd[:, 1] = posts[:, 0]
        
        cmd_plot(cmd, pcmd, ['DR10', 'XD Post.'], os.environ['xdplots'] + 'foo.png')

"""
       single = np.vstack((Xsingle[:, 6] + Xsingle[:, 0],
                            Xsingle[:, 7] + Xsingle[:, 0],
                            Xsingle[:, 0],
                            Xsingle[:, 8] + Xsingle[:, 0],
                            Xsingle[:, 9] + Xsingle[:, 0])).T

        coadd = np.vstack((Xcoadd[:, 6] + Xcoadd[:, 0],
                            Xcoadd[:, 7] + Xcoadd[:, 0],
                            Xcoadd[:, 0],
                            Xcoadd[:, 8] + Xcoadd[:, 0],
                            Xcoadd[:, 9] + Xcoadd[:, 0])).T

        post = np.vstack((posts[:, 6] + posts[:, 0],
                          posts[:, 7] + posts[:, 0],
                          posts[:, 0],
                          posts[:, 8] + posts[:, 0],
                          posts[:, 9] + posts[:, 0])).T
"""
