# Author - Ross Fadely
#
import numpy as np
np.seterr(divide='ignore')
import ctypes as ct
import pyfits as pf

from scipy.optimize import fmin_l_bfgs_b
from utils import *


class HBsep(object):
    """
    Hierarchical Bayesian classification of astronomical 
    objects, using photometry.
    """
    def __init__(self, class_labels, Nzs, z_maxs=None, z_min=0.):

        self.Nzs = Nzs
        self.z_min = z_min
        
        if z_maxs is None:
            self.z_maxs = np.zeros(len(Nzs))
        else:
            self.z_maxs = z_maxs
            
        self.class_labels = class_labels
        self.Nclasses = len(class_labels)

        self.model_mags = {}
        self.model_fluxes = {}

    def get_filter_norm(self, filter_list_path):
        """
        Return the normalizing flux in AB for given
        filter list.
        """
        f = open(filter_list_path)
        self.Nfilter = len(f.readlines())
        f.close
        self._make_models(None, filter_list_path, 1, 0., 0.,
                          np.zeros((2, 2)).astype(np.float64), True)

    def _make_models(self, sed_list_path, filter_list_path, Nz, zmin, zmax,
                     models, filter_only=False):
        """
        Prepare for and call model maker.
        """
        # load funtion
        model_maker = ct.CDLL('./_model_maker.so')

        # initialize filter normalizations
        self.filter_norms = np.zeros(self.Nfilter).astype(np.float64)

        # pointers for filters
        filt_norm_p = self.filter_norms.ctypes.data_as(ct.POINTER(ct.c_double))
        filt_list_p = ct.c_char_p(filter_list_path)

        # ctypes foo
        Nz = ct.c_long(np.int64(Nz))
        zmin = ct.c_double(np.float64(zmin))
        zmax = ct.c_double(np.float64(zmax))
        models_p = ctype_2D_double_pointer(models)

        # filters only, or models too?
        if filter_only:
            sed_list_p = filt_list_p
            filter_only = ct.c_long(np.int64(1))
        else:
            sed_list_p = ct.c_char_p(sed_list_path)
            filter_only = ct.c_long(np.int64(0))

        model_maker.model_maker(filt_list_p, sed_list_p, filter_only, Nz, zmin,
                                zmax, filt_norm_p, models_p)

    def create_models(self, filter_list_path, list_of_sed_list_paths,
                      normalize_models=True):
        """
        For each class, produce models over redshifts
        """

        # specify Nfilter
        f = open(filter_list_path)
        self.Nfilter = len(f.readlines())
        f.close

        for i in range(self.Nclasses):
            key = self.class_labels[i]
            sed_list_path = list_of_sed_list_paths[i]

            # get number of seds in class
            f = open(sed_list_path)
            Nseds = len(f.readlines())
            f.close

            # init
            self.model_fluxes[key] = \
                np.zeros((self.Nzs[i] * Nseds,
                          self.Nfilter)).astype(np.float64)

            # make the models
            self._make_models(sed_list_path, filter_list_path, self.Nzs[i],
                              self.z_min, self.z_maxs[i],
                              self.model_fluxes[key])

            # normalize models
            if normalize_models:
                self.model_fluxes[key] /= np.mean(self.model_fluxes[key],
                                                  axis=1)[:, None]

            # compute magnitudes
            self.model_mags[key] = -2.5 * np.log10(self.model_fluxes[key] /
                                                   self.filter_norms[None, :])

            print '\nCreated '+key+' models'

    def read_and_process_data(self, data,
                              missing_mags=None,
                              limiting_mags=None,
                              limiting_sigmas=None,
                              normalize_flux=True,
                              inflation_factor=1e6):
        """
        Read in data location, process.
        """
        self.missing_mags = missing_mags
        self.limiting_mags = limiting_mags
        self.limiting_sigmas = limiting_sigmas

        # read and process data
        data = self.read_data(data)
        self.Ndata = data.shape[0]
        self.mags = data[:, :self.Nfilter]
        self.mag_errors = data[:, self.Nfilter:]
        self.process_data(inflation_factor, normalize_flux)

    def read_data(self, data):
        """
        Read in data file or array, do some moderate
        error checking.
        """
        if isinstance(data, np.ndarray):
            pass
        elif isinstance(data, str):
            try:
                # hdu number might be an issue
                f = pf.open(data)
                tdata = f[1].data
                f.close()
                data = np.zeros((len(tdata.field(0)), len(tdata[0])))

                for m in range(data.shape[1]):
                    data[:, m] = tdata.field(m)
            except:
                try:
                    data = np.loadtxt(data)
                except:
                    print '\nData file not read.'
                    assert False

        self.shape_check(data, 'Data Error:')
        return data

    def shape_check(self, data, item):
        """
        Simple shape checks
        """
        assert len(data.shape) == 2, '\n\n' + item + \
            ' Must be a 2D numpy array'
        assert data.shape[0] > data.shape[1], '\n\n' + item + \
            ' Must have more rows than columns'
        assert np.mod(data.shape[1], 2) == 0, '\n\n' + item + \
            ' Ncolumn is not even, what gives?'

    def process_data(self, inflation_factor, normalize_flux):
        """
        Turn data into flux, flux err
        """
        # must calculate normalization first
        try:
            print '\nFilter AB normalizations: ', self.filter_norms
            print
        except:
            assert False, 'Must run make_model or get_filter_norm'

        # normal flux calculation
        self.fluxes = 10.0 ** (-0.4 * self.mags) * self.filter_norms[None, :]
        self.flux_errors = 10.0 ** (-0.4 * self.mags) * \
            self.filter_norms[None, :] * np.log(10.0) * 0.4 * self.mag_errors

        # account for missing data
        if self.missing_mags is not None:
            Nind, find = np.where(self.mags == self.missing_mags)
            ind = np.where(self.fluxes > 0)[0]
            self.fluxes[Nind,find] = np.median(self.fluxes[ind])
            self.flux_errors[Nind,find] = self.fluxes[Nind,find] * \
                inflation_factor

        # account for anything fainter than limiting magnitudes
        if self.limiting_mags is not None:
            assert self.limiting_sigmas is not None, \
                'Must specify Nsigma for limiting mags'
            for i in range(self.Nfilter):
                ind = np.where(self.mags[:,i] >= self.limiting_mags[i])[0]
                self.fluxes[ind, i] = 0.0
                self.flux_errors[ind, i] = 10.0**(-0.4*self.limiting_mags[i]) \
                    * self.limiting_sigmas[i]

        # normalize if desired
        if normalize_flux:
            ind = np.where(self.fluxes != 0.0)[0]
            self.flux_norm = np.mean(self.fluxes[ind])
            self.fluxes /= self.flux_norm
            self.flux_errors /= self.flux_norm

    def fit_models(self):
        """
        Produce chi2 and coefficients of model fits to data.
        """
        fitter = ct.CDLL('./_fit_models.so')

        self.chi2s = {}
        self.coeffs = {}
        self.coefferrs = {}
        for i in range(self.Nclasses):

            key = self.class_labels[i]
            Nmodel = self.model_fluxes[key].shape[0]
            self.chi2s[key] = np.zeros((self.Ndata, Nmodel))
            self.coeffs[key] = np.zeros((self.Ndata, Nmodel))
            self.coefferrs[key] = np.zeros((self.Ndata, Nmodel))

            # ctypes prep
            Ndata = ct.c_long(self.Ndata)
            Nmodel = ct.c_long(Nmodel)
            Nfilter = ct.c_long(self.Nfilter)
            datap = ctype_2D_double_pointer(self.fluxes)
            dataerrp = ctype_2D_double_pointer(self.flux_errors)
            modelsp = ctype_2D_double_pointer(self.model_fluxes[key])
            chi2sp = ctype_2D_double_pointer(self.chi2s[key])
            coeffsp = ctype_2D_double_pointer(self.coeffs[key])
            coefferrsp = ctype_2D_double_pointer(self.coefferrs[key])

            fitter.fit_models(Ndata, Nmodel, Nfilter, modelsp, datap, dataerrp,
                              coeffsp, coefferrsp, chi2sp)

    def coefficient_marginalization(self, Nstep_factor=2, Nsigma=5,
                                    delta_chi2_cut=32., floor=1e-100):
        """
        Marginalize over the fit coefficients.
        """
        Ndata = ct.c_long(self.Ndata)
        Nstep = ct.c_long((Nsigma * Nstep_factor) * 2 + 1)
        Nsigma = ct.c_double(Nsigma)
        det_flux_errors = np.prod(self.flux_errors, axis=1)
        det_flux_errorsp = det_flux_errors.ctypes.data_as(
            ct.POINTER(ct.c_double))
        delta_chi2_cut = ct.c_double(delta_chi2_cut)

        marg = ct.CDLL('./_coeff_marginalization.so')
        self.calc_coeff_priors()

        self.ignored = np.zeros(self.Ndata)
        self.bad_fit_flags = {}
        self.coeff_marg_like = {}
        for i in range(self.Nclasses):

            key = self.class_labels[i]
            Nmodel = self.model_fluxes[key].shape[0]
            minchi2 = np.min(self.chi2s[key], axis=1)
            self.coeff_marg_like[key] = np.zeros((self.Ndata, Nmodel))

            # ctypes prep
            Nmodel = ct.c_long(Nmodel)
            minchi2p = minchi2.ctypes.data_as(ct.POINTER(ct.c_double))
            prior_varsp = self.coeff_prior_vars[key].ctypes.data_as(
                ct.POINTER(ct.c_double))
            prior_meansp = self.coeff_prior_means[key].ctypes.data_as(
                ct.POINTER(ct.c_double))
            chi2sp = ctype_2D_double_pointer(self.chi2s[key])
            coeffsp = ctype_2D_double_pointer(self.coeffs[key])
            coefferrsp = ctype_2D_double_pointer(self.coefferrs[key])
            marglikep = ctype_2D_double_pointer(self.coeff_marg_like[key])

            marg.coeff_marginalization(Nstep, Ndata, Nmodel, Nsigma, minchi2p,
                                       delta_chi2_cut, coeffsp, coefferrsp,
                                       prior_meansp, prior_varsp,
                                       det_flux_errorsp, chi2sp, marglikep)

            # Flag bad fits
            self.bad_fit_flags[key] = np.zeros(self.Ndata)
            ind = np.where(self.coeff_marg_like[key].sum(axis=1) < floor)[0]
            self.bad_fit_flags[key][ind] = 1
            self.ignored[ind] += 1

        # construct index for likelihood estimation
        ind = np.where(self.ignored < self.Nclasses)[0]
        self.use = ind
        ind = np.where(self.ignored == self.Nclasses)[0]
        self.ignore = ind

    def calc_coeff_priors(self):
        """
        Compute prior parms for coefficient fits.
        """
        self.coeff_prior_vars = {}
        self.coeff_prior_means = {}
        for i in range(self.Nclasses):

            key = self.class_labels[i]
            coeffs = self.coeffs[key]
            weights = 1./self.coefferrs[key]

            mean = np.sum(weights * np.log(coeffs), axis=0) / \
                np.sum(weights, axis=0)
            var = np.var(coeffs, axis=0)
            self.coeff_prior_means[key] = mean
            self.coeff_prior_vars[key] = var

    def apply_and_marg_redshift_prior(self, method):
        """
        Apply redshift prior and margninalize over redshift
        """
        self.zc_marg_like = {}
        for i in range(self.Nclasses):
            key = self.class_labels[i]

            Nz = self.Nzs[i]
            Nmodel = self.model_fluxes[key].shape[0]
            Ntemplate = Nmodel/Nz

            if Nz == 1:
                self.zc_marg_like[key] = self.coeff_marg_like[key]
                continue
            zgrid = np.linspace(0, self.z_maxs[i], Nz)

            self.zc_marg_like[key] = np.zeros((self.Ndata, Ntemplate))
            if method == 1:
                # shape = Nmodel,Nz
                z_medians = np.array([np.ones(Nz) * zmed
                                      for zmed in self.z_medians[key]])
                z_prior = zgrid**self.z_pow[key] * \
                    np.exp(-(zgrid/z_medians)**self.z_pow[key])
                z_prior /= z_prior.sum(axis=1)[:, None]
                prior_weighted_like = self.coeff_marg_like[key] * \
                    z_prior.ravel()[None, :]
                for j in range(Ntemplate):
                    self.zc_marg_like[key][:, j] = \
                        np.sum(prior_weighted_like[:, j*Nz:j*Nz+Nz], axis=1)

    def assign_hyperparms(self, hyperparms, method):
        """
        Assign hyperparameters from a flattened list (from optimizer)
        """
        count = 0

        # assign template weights
        self.template_weights = {}
        for i in range(self.Nclasses):
            key = self.class_labels[i]
            Ntemplate = np.int(self.model_fluxes[key].shape[0] / self.Nzs[i])
            self.template_weights[key] = \
                np.exp(np.array(hyperparms[count:count+Ntemplate]))
            self.template_weights[key] /= self.template_weights[key].sum()
            count += Ntemplate

        # prior parms
        if method == 1:
            self.z_medians = {}
            for i in range(self.Nclasses):
                if self.Nzs[i] == 1:
                    continue
                key = self.class_labels[i]
                self.z_medians[key] = \
                    np.array(hyperparms[count:count+Ntemplate])
                count += Ntemplate
            self.z_pow = {}
            for i in range(self.Nclasses):
                if self.Nzs[i] == 1:
                    continue
                self.z_pow[key] = hyperparms[count:count+1]
                count += 1

        self.class_weights = np.exp(np.array(hyperparms[-self.Nclasses:]))
        self.class_weights /= self.class_weights.sum()

    def calc_neg_lnlike(self, method, floor=1e-100):
        """
        Calculate marginalized likelihoods.
        """
        self.tzc_marg_like = {}
        self.marg_like = np.zeros(self.Ndata)
        ind = self.use
        for i in range(self.Nclasses):
            key = self.class_labels[i]
            self.tzc_marg_like[key] = np.zeros(self.Ndata)
            self.tzc_marg_like[key][ind] = \
                np.sum(self.zc_marg_like[key][ind] *
                       self.template_weights[key][None, :], axis=1)

            self.marg_like += np.maximum(self.tzc_marg_like[key] *
                                         self.class_weights[i], floor)

        self.neg_log_likelihood = -1.0 * np.sum(np.log(self.marg_like[ind]))

    def call_neg_lnlike(self, hyperparms, method=1):
        """
        Give this to optimizer to call.
        """
        weights = np.exp(np.array(hyperparms[-self.Nclasses:]))
        if np.Inf in weights:
            return np.Inf
        self.assign_hyperparms(hyperparms, method)
        self.apply_and_marg_redshift_prior(method)
        self.calc_neg_lnlike(method)
        return self.neg_log_likelihood

    def init_hyperparms(self, z_median, z_pow):
        """
        Initialize flattened list of hyperparameters.
        """
        # initialize parameters, barf
        p0 = np.array([])
        for i in range(self.Nclasses):
            key = self.class_labels[i]
            Ntemplate = self.model_fluxes[key].shape[0] / self.Nzs[i]
            p0 = np.append(p0, np.ones(Ntemplate) * np.log(1./Ntemplate))
        for i in range(self.Nclasses):
            if self.Nzs[i] != 1:
                p0 = np.append(p0, np.ones(Ntemplate) * z_median[i])
        for i in range(self.Nclasses):
            if self.Nzs[i] != 1:
                p0 = np.append(p0, z_pow[i])
        for i in range(self.Nclasses):
            p0 = np.append(p0, np.log(1./self.Nclasses))

        return p0

    def init_hyperparm_bounds(self):
        """
        Make bounds for fmin_l_bfgs_b
        """
        bounds = []
        for i in range(self.Nclasses):
            key = self.class_labels[i]
            Ntemplate = self.model_fluxes[key].shape[0] / self.Nzs[i]
            bounds.extend([(-1.*np.Inf, np.Inf) for j in range(Ntemplate)])
        for i in range(self.Nclasses):
            if self.Nzs[i] != 1:
                bounds.extend([(0.05, self.z_maxs[i]) for j in
                               range(Ntemplate)])
        for i in range(self.Nclasses):
            if self.Nzs[i] != 1:
                bounds.extend([(0., 4)])
        for i in range(self.Nclasses):
            bounds.extend([(-1.*np.Inf, np.Inf)])

        return bounds

    def optimize(self, z_median=None, z_pow=None, init_p0=None,
                 eps=1.e-1, factr=1.e7, maxfun=1):
        """
        Optimize using scipy's fmin_l_bfgs_b
        """
        if init_p0 is not None:
            p0 = init_p0
        else:
            p0 = self.init_hyperparms(z_median, z_pow)

        bounds = self.init_hyperparm_bounds()

        self.init_nll = self.call_neg_lnlike(p0)
        result = fmin_l_bfgs_b(self.call_neg_lnlike, p0, approx_grad=1,
                               bounds=bounds, epsilon=eps, factr=factr,
                               maxfun=maxfun, iprint=2)

        return self.call_neg_lnlike(result[0]), result[0]

    def get_relative_likelihoods(self):
        """
        Return array of relative likelihoods
        """
        relative_likelihoods = np.zeros((self.Ndata, self.Nclasses))
        for i in range(self.Nclasses):
            key = self.class_labels[i]
            relative_likelihoods[:, i] = self.tzc_marg_like[key] * \
                self.class_weights[i]

        relative_likelihoods /= relative_likelihoods.sum(axis=1)[:, None]
        return relative_likelihoods

    def write_fits_table(self, filename, data, labels):
        """
        Write a FITS data table to given `filename` usings `labels` for 
        column names.
        """
        Nlabels = len(labels)
        if Nlabels > 1:
            cols = pf.ColDefs([pf.Column(name=labels[i], format='E',
                                         array=data[:, i]) 
                               for i in range(Nlabels)])
        else:
            cols = pf.ColDefs([pf.Column(name=labels[0], format='E',
                                         array=data[0,:])])

        tbhdu = pf.new_table(cols)
        hdu = pf.PrimaryHDU(np.arange(10))

        tblist = pf.HDUList([hdu, tbhdu])
        tblist.writeto(filename, clobber=True)

    def write_relative_likelihoods(self, filename):
        """
        Write relative marg. likelihoods to FITS file.
        """
        likes = self.get_relative_likelihoods()
        flags = np.ones(likes.shape[0])

        ind = np.in1d(self.ignore, np.arange(likes.shape[0]))
        flags[ind] = 0.0

        out = np.zeros((likes.shape[0], likes.shape[1]+1))
        out[:, :-1] = likes
        out[:, -1] = flags
        labels = [self.class_labels[i] for i in range(self.Nclasses)]
        labels.append('used')
        self.write_fits_table(filename, out, labels)

    def write_array(self, filename, array, label):
        """
        Write a 1D array to FITS table.
        """
        self.write_fits_table(filename, np.atleast_2d(array), label)

    def write_minchi2(self, filename):
        """
        Write the min. chi^2 for each class to FITS table.
        """
        labels = []
        minchi2s = np.zeros((self.Ndata, self.Nclasses))
        for i in range(self.Nclasses):
            key = self.class_labels[i]
            labels.append(key+' minchi2')
            minchi2s[:, i] = self.chi2s[key].min(axis=1)

        self.write_fits_table(filename, minchi2s, labels)
