#
# Plotting routines, mostly for exploration.  Paper figures are generated with `fig_x`
#
# Author - Ross Fadely
#
from matplotlib import use; use('Agg')
from utils import fetch_mixed_epoch

import numpy as np
import matplotlib.pyplot as pl

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

if __name__ == '__main__':
    #one_epoch_class_check()
    plot_misclass()
