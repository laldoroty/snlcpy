import pycmpfit
import pandas as pd
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import shape
import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d
import os
import sys
# sys.path.insert(1, '/home/astrolab/pycmpfit/build/lib.linux-x86_64-3.6')
# from kapteyn import kmpfit

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


def fit_pcparams(*data, fpca_f='vague', init_guess=None,
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
            preidx = np.sum(model_x < penalty_increase)
            pre_diff = model_y[1:preidx] - model_y[0:preidx-1]
            pre_violate = np.sum(np.where(pre_diff < 0, 0, pre_diff))

            # thrid add the constraints of monotonically decrease after penalty_decrease
            afteridx = np.sum(model_x < penalty_decrease)
            after_diff = model_y[afteridx+1:] - model_y[afteridx:-1]
            after_violate = np.sum(np.where(after_diff > 0, 0, after_diff))

            resid += max_violate*10 + after_violate * \
                10 + abs(pre_violate*10)    #

        user_dict = {"deviates": resid}
        return user_dict

    if init_guess is None:
        init_guess = [date[np.argmin(mag)], 17, 1, 0, 0, 0]
    init_guess = np.array(init_guess)

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
    :params mpfit_result
    :params fpca_dir: directory where fpca templates are at
    :params return_func: if true, return function, if not, return grids

    return function or grids depending on return_func
    '''
    if fpca_f in ['B', 'V', 'R', 'I']:
        date = np.arange(-10, 50, 0.1)
    else:
        date = np.arange(-10, 40, 0.1)

    PCvecs = get_pctemplates(fpca_f)
    coeffs = np.array([fit_result['params'][i] for i in range(2, 6)])
    coeffs = np.insert(coeffs, 0, 1)
    print(coeffs)
    date += fit_result['params'][0]
    print(PCvecs)
    print(len(PCvecs))
    LC = fit_result['params'][1] + np.sum(np.dot(PCvecs, coeffs))
    lc_error = None
    # all_the_jacobians = np.dstack((np.zeros(len(date)),
    #     np.zeros(len(date)),
    #     np.full(len(date),PCvecs[0]),
    #     np.full(len(date),PCvecs[1]),
    #     np.full(len(date),PCvecs[2]),
    #     np.full(len(date),PCvecs[3])))

    # lc_error = np.sqrt(np.matmul(np.matmul(all_the_jacobians,
    #     fit_result['mpfit_result'].covar)),np.transpose(all_the_jacobians))

    if return_func == False:
        return date, LC, lc_error
    else:
        return interp1d(date, LC), interp1d(date, lc_error)

###############################################################################
# # test_dict = {'date': [1, 2, 3], 'mag': [3, 4, 5], 'emag': [4, 6, 3]}
# # fit_pcparams(test_dict)
# test_tab = Table([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
#                  names=('date', 'mag', 'emag'))
# print('success')
# fit_pcparams(test_tab)
# print('success #2')


fpca_f = 'B'

basis = get_pctemplates(fpca_f)
# print(basis[0](-20))


def get_modelval(date, theta, fpca_f):
    '''
    Calculates the linear combination of coefficients and PC vectors at every point in the grid.
    Given a phase and parameters, returns the corresponding fitted magnitude.
    :params phase:
    :params theta:
    '''
    basis = get_pctemplates(fpca_f)
    if not isinstance(date, float):
        date = np.array(date)
    tmax, mmax, a1, a2, a3, a4 = theta
    coeffs = [1, a1, a2, a3, a4]
    mphase = date-tmax
    y = mmax + (np.dot(coeffs,
                       np.array([fbasis(mphase) for fbasis in basis])))

    return y


theta = [0, 0, 10, 0, 0, 0]
x = np.arange(-10, 60, 1)
# print(get_modelval([0, 1], theta))
# plt.scatter(x, get_modelval(x, theta))
# plt.show()

a = pd.read_csv(os.path.join(abspath, '02boPriv.csv'))
a = a[a.Passband == 'V(kait3)']
# print(a)
data = {'date': np.array(a['MJD_OBS']) - a['MJD_OBS'].tolist()[np.argmin(np.array(a['MAG']))], 'mag': np.array(
    a['MAG']), 'emag': np.array(a['eMAG'])}

res = fit_pcparams(data, fpca_f='vague', init_guess=None,
                   components=2, penalty=True, penalty_increase=None,
                   penalty_decrease=None, boundary=None)

print('parameters', res['params'])
# print(make_fittedlc(fpca_f, res, fpca_dir='', return_func=True))
print(get_modelval(data['date'], res['params'], 'vague'))
plt.scatter(data['date'], get_modelval(
    data['date'], res['params'], 'vague'), label='model')
plt.scatter(data['date'], data['mag'], label='data')
plt.legend()
plt.show()
