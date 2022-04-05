from unittest import result
import pandas as pd
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import shape
import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d
import os
import sys
import pycmpfit

abspath = os.path.dirname(__file__)
FPCA_dir = os.path.join(abspath, 'LCfPCA_He')


# def phase_limit(fpca_f):
#     '''
#     :params fpca_f: Template filter that you'd like to use. B, V, R, I filters available, Otherwise, use band vague
#     '''
#     if fpca_f in ['B', 'V', 'R', 'I']:
#         phaselim = [-10, 50]
#     else:
#         phaselim = [-10, 40]

#     return phaselim


def get_pctemplates(fpca_f):
    '''
    returns corresponding fpca templates
    :params fpca_f: Template filter that you'd like to use. B, V, R, I filters available, Otherwise, use band vague template.
    return: list of interpolated template functions [mean, pc1, pc2, pc3, pc4] -- note that returned interpolated functions\
                    are in -10 to 40/50 for time axis (40 for band vague, 50 for BVRI)
    '''
    if fpca_f in ['B', 'V', 'R', 'I']:
        fname = 'bandSpecific_%s.txt' % fpca_f
    else:
        fname = 'bandVague.txt'

    fpca_file = os.path.join(FPCA_dir, fname)  # FIXME
    t_pc = Table.read(fpca_file, format='ascii')
    phase, mean, pc1, pc2, pc3, pc4 = t_pc['phase'], t_pc['mean'], \
        t_pc['FPC1'], t_pc['FPC2'], t_pc['FPC3'], t_pc['FPC4']

    fpc0 = interp1d(phase, mean, fill_value=np.nan, bounds_error=False)
    fpc1 = interp1d(phase, pc1, fill_value=np.nan, bounds_error=False)
    fpc2 = interp1d(phase, pc2, fill_value=np.nan, bounds_error=False)
    fpc3 = interp1d(phase, pc3, fill_value=np.nan, bounds_error=False)
    fpc4 = interp1d(phase, pc4, fill_value=np.nan, bounds_error=False)

    return [fpc0, fpc1, fpc2, fpc3, fpc4]


def fit_pcparams(data, fpca_f='vague', init_guess=None,
                 components=2, penalty=True, penalty_increase=None,
                 penalty_decrease=None, boundary=None):
    '''
    :params data: list of dates, OR 3*n array or n*3 array when mag and magerr are None
    :params fpca_f: fpca filter that you want to use to fit the lc. options:\
                    vague, B, V, R, I
    :params init_guess: list, initial guess of parameters, None for default
    :params components: number of fPCA components you want to use
    :params penalty:
    :params penalty_increase: latest epoch before which you'd like the lc to be monotonically increasing
    :params penalty_decrease: earliest epoch after which you'd like the lc  to be monotonically decreasing
    :params boundary: lists of shape components*2 [[low, high], [low, high]], put constraints on pc parameters

    return fpca_f, mpfit_result
    mpfit_result: as long as it contains attributes params, covar, chi2
    '''
    # Read in data if given separate arrays or lists
    if len(data) == 3:
        date, mag, emag = [np.array(ii) for ii in data]
    if len(data) == 2:
        date, mag = [np.array(ii) for ii in data]
        # if no magnitude error is provided, then use constant 1 for all epochs.
        # This makes the denominator in the merit function = 1
        emag = np.array([1]*len(date))

    # Read in data if given a dictionary
    if len(data) == 1 and isinstance(data[0], dict):
        data = data[0]
        date, mag = map(lambda x: np.array(data[x]), ['date', 'mag'])
        # I think this should be data.keys() in the next line:
        if 'emag' not in data:
            emag = np.array([1]*len(date))
        else:
            emag = np.array(data['emag'])

    # Read in data if given an astropy table
    if len(data) == 1 and isinstance(data[0], Table):
        data = data[0]
        date, mag = np.array(data['date']), np.array(data['mag'])
        if 'emag' not in data.colnames:
            emag = np.array([1]*len(date))
        else:
            emag = np.array(data['emag'])

    def get_modelval(date, theta, fpca_f):
        '''
        Calculates the linear combination of coefficients and PC vectors at every point in the grid.
        Given a phase and parameters, returns the corresponding fitted magnitude.
        :params date: float or list/array
        :params theta: list of length 6. fitting parameters 

        return: float if input is float, array if input is array [tmax, mmax, a1, a2, a3, a4]
        '''
        basis = get_pctemplates(fpca_f)
        if not isinstance(date, float):
            date = np.array(date)
        tmax, mmax, a1, a2, a3, a4 = theta
        coeffs = [1, a1, a2, a3, a4]
        mphase = date-tmax
        y = mmax + (np.dot((coeffs),
                           np.array([fbasis(mphase) for fbasis in basis])))

        return y
    
    def meritfunc(m, n, theta, input_data):
        # , penalty=penalty, \
        #       penalty_increase=penalty_increase, penalty_decrease=penalty_decrease):
        # FIXME not sure if mpfit allows users to add more arguments in userfunc?
        '''
        Chisq to be minimized.
        :m number of samples (len(data)):
        :n number of parameters (len(theta)):
        '''
        x, y, ey = map(lambda xx: np.array(
            input_data[xx]), ['date', 'mag', 'emag'])
        fpca_f, penalty_increase, penalty_decrease = map(
            lambda x: input_data[x], ['fpca_f', 'penalty_increase', 'penalty_decrease'])
        y_model = get_modelval(x, theta, fpca_f)
        resid = (y - y_model)/ey
        resid = np.where(np.isnan(resid), 0, resid)
        # resid_dict
        if penalty:
            if penalty_increase is None:
                penalty_increase = -1  # -0.03? FIXME
            if penalty_decrease is None:
                penalty_decrease = 35

            model_x = np.arange(-10, 50, 0.01) + theta[0]
            # print(model_x)
            model_y = get_modelval(model_x, theta, fpca_f)
            idx = ~np.isnan(model_y)
            model_x = model_x[idx]
            model_y = model_y[idx]

            # first let's make sure the fitted Tmax is the real Tmax.
            # --> this is added bc when coefficients of higher PC templates
            # are too large, the acutal maximum on the fitted line maybe
            # at different epoch that the Tmax we fitted.

            # True maxmimum magnitude at fitted Tmax
            max_mag = get_modelval(theta[0], theta, fpca_f)
            # or max_mag = get_modelval(-0.03?, theta) FIXME

            # whenever fitted lc is brighter than maxmimum date's magnitude, add that difference to residual.
            max_violate = np.sum(abs(model_y-max_mag)[model_y < max_mag])

            # second add the constraints of monotonically increase before penalty_increase
            # note preidx here is a number
            preidx = np.sum(model_x < penalty_increase+theta[0])
            pre_diff = model_y[1:preidx] - model_y[0:preidx-1]
            pre_violate = np.sum(np.where(pre_diff < 0, 0, pre_diff))

            # thrid add the constraints of monotonically decrease after penalty_decrease
            afteridx = np.sum(model_x < penalty_decrease)
            after_diff = model_y[afteridx+1:] - model_y[afteridx:-1]
            after_violate = np.sum(np.where(after_diff > 0, 0, after_diff))

            resid += max_violate*10 + abs(after_violate) * \
                10 + abs(pre_violate*10)    #

        user_dict = {"deviates": resid}
        return user_dict

    if init_guess is None:

        init_guess = [date[np.argmin(mag)], 17, 1, 0, 0, 0]
    init_guess = np.array(init_guess)
    # print(init_guess)
    m, n = len(date), len(init_guess)
    input_data = {'date': date, 'mag': mag, 'emag': emag,
                  'fpca_f': fpca_f, 'penalty_increase': penalty_increase, 'penalty_decrease': penalty_decrease}

    py_mp_par = list(pycmpfit.MpPar() for i in range(n))
    py_mp_par[0].limited = [1, 1]
    py_mp_par[0].limits = [init_guess[0]-10, init_guess[0]+10]

    for ii in range(4-components):
        py_mp_par[-(ii+1)].fixed = True
    if boundary:
        for ii in range(components):
            py_mp_par[2+ii].limited = [1, 1]
            py_mp_par[2+ii].limits = boundary[ii]
    # print(input_data)
    # print(init_guess)
    fit = pycmpfit.Mpfit(meritfunc, m, init_guess,
                         private_data=input_data, py_mp_par=py_mp_par)
    fit.mpfit()  # NOTE now init_guess has been updated

    fit_result = {'mpfit_result': fit.result, 'params': init_guess, }

    # calculate true reduced chi square
    # resid = get_modelval()
    return fit_result


def make_fittedlc(fpca_f, fit_result, fpca_dir='', return_func=True):
    '''
    :params fpca_f: fpca filter, options: vague, B, V, R, I
    :params mpfit_result: output from fit_pcparams()
    :params fpca_dir: directory where fpca templates are at
    :params return_func: if true, return function, if not, return grids

    return function or grids depending on return_func
    '''
    if fpca_f in ['B', 'V', 'R', 'I']:
        date = np.arange(-10, 50, 0.1)
    else:
        date = np.arange(-10, 40, 0.1)
    
    PCvecs = np.array(get_pctemplates(fpca_f))
    PCvecs_discrete = []
    for vec in PCvecs:
        PCvecs_discrete.append([vec(t) for t in date])
    PCvecs_discrete = np.array(PCvecs_discrete)

    coeffs = np.array([fit_result['params'][i] for i in range(2, 6)])
    coeffs = np.insert(coeffs, 0, 1)
    date += fit_result['params'][0]
    LC = fit_result['params'][1] + np.dot(np.transpose(PCvecs_discrete), coeffs)
    all_the_jacobians = np.column_stack((np.zeros(len(date)),
        np.ones(len(date)),
        PCvecs_discrete[0],
        PCvecs_discrete[1],
        PCvecs_discrete[2],
        PCvecs_discrete[3]))

    lc_error = []
    # I don't want this to be a loop, but it works, so make this a vector operation later. 
    for jac in all_the_jacobians:
        j = np.sqrt(np.matmul(jac,np.matmul(fit_result['mpfit_result'].covar,np.transpose(jac))))
        lc_error.append(j)
    lc_error = np.array(lc_error)

    if return_func == False:
        return date, LC, lc_error
    else:
        return interp1d(date, LC, bounds_error=False), interp1d(date, lc_error, bounds_error=False)

def plot_lc(data,fit_result,fpca_f, fpca_dir='', input_func=True):
    '''
    Shortcut to plot light curves. 
    :params input_func: True if input is continuous, i.e., scipy.interpolate.interp1d(). False if arrays are provided.
    :params data
    :params fit: Result from make_fittedlc().
    '''

    fig, ax =plt.subplots(figsize=(10,8))
    fit = make_fittedlc(fpca_f, fit_result,fpca_dir=fpca_dir, return_func=True)

    x = np.arange(-10.0,50.0,0.1)+fit_result['params'][0]

    # Read in data if given separate arrays or lists
    if len(data) == 3:
        date, mag, emag = [np.array(ii) for ii in data]
    if len(data) == 2:
        date, mag = [np.array(ii) for ii in data]
        # if no magnitude error is provided, then use constant 1 for all epochs.
        # This makes the denominator in the merit function = 1
        emag = np.array([1]*len(date))

    # Read in data if given a dictionary
    if len(data) == 1 and isinstance(data[0], dict):
        data = data[0]
        date, mag = map(lambda x: np.array(data[x]), ['date', 'mag'])
        if 'emag' not in data:
            emag = np.array([1]*len(date))
        else:
            emag = np.array(data['emag'])

    # Read in data if given an astropy table
    if len(data) == 1 and isinstance(data[0], Table):
        data = data[0]
        date, mag = np.array(data['date']), np.array(data['mag'])
        if 'emag' not in data.colnames:
            emag = np.array([1]*len(date))
        else:
            emag = np.array(data['emag'])

    ax.scatter(date,mag,s=20,edgecolor='k')
    ax.errorbar(date,mag,yerr=emag,color='C0',markersize=10,ls='none')

    if input_func: 
        ax.plot(x, fit[0](x), color='C0', label='Fit')
        ax.fill_between(x, fit[0](x) - fit[1](x), fit[0](x) + fit[1](x), color='C0', alpha=0.5)
    else:
        ax.plot(fit[0], fit[1], color='C0', label='Fit')
        ax.fill_between(fit[0], fit[1] - fit[0], fit[1] + fit[0], color='C0', alpha=0.5)

    ax.invert_yaxis()
    plt.legend()
    plt.show()
