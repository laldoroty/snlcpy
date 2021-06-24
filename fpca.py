import pandas as pd
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import shape
import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d
import os
# from kapteyn import kmpfit
FPCA_dir = './LCfPCA_He'


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

    fpca_file = os.path.join(FPCA_dir, fname)
    t_pc = Table.read(fpca_file, format='ascii')
    phase, mean, pc1, pc2, pc3, pc4 = t_pc['phase'], t_pc['mean'], \
        t_pc['FPC1'], t_pc['FPC2'], t_pc['FPC3'], t_pc['FPC4']

    fpc0 = interp1d(phase, mean, fill_value=np.nan, bounds_error=False)
    fpc1 = interp1d(phase, pc1, fill_value=np.nan, bounds_error=False)
    fpc2 = interp1d(phase, pc2, fill_value=np.nan, bounds_error=False)
    fpc3 = interp1d(phase, pc3, fill_value=np.nan, bounds_error=False)
    fpc4 = interp1d(phase, pc4, fill_value=np.nan, bounds_error=False)

    return [fpc0, fpc1, fpc2, fpc3, fpc4]


def fit_pcparams(*data, mag=None, magerr=None, fpca_f='vague', init_guess=None,
                 components=2, fpca_dir='', penalty=True, penalty_increase=None,
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

    # Initial guess
    init_maxdate = date[np.argmin(mag)]

    basis = get_pctemplates(fpca_f)

    def get_modelval(date, theta):
        '''
        Calculates the linear combination of coefficients and PC vectors at every point in the grid.
        Given a phase and parameters, returns the corresponding fitted magnitude.
        :params date: float or list/array
        :params theta: list of length 6. fitting parameters 

        return: float if input is float, array if input is array [tmax, mmax, a1, a2, a3, a4]
        '''
        if not isinstance(date, float):
            date = np.array(date)
        tmax, mmax, a1, a2, a3, a4 = theta
        coeffs = [1, a1, a2, a3, a4]
        mphase = date-tmax
        y = mmax + (np.dot((coeffs),
                           np.array([fbasis(mphase) for fbasis in basis])))

        return y

    def meritfunc(theta, input_data):
        '''
        Chisq to be minimized.
        :m number of samples (len(data)):
        :n number of parameters (len(theta)):
        '''
        x, y, ey = map(lambda xx: np.array(
            input_data[xx]), ['date', 'mag', 'emag'])
        y_model = get_modelval(x, theta)
        resid = (y - y_model)/ey
        resid = np.where(np.isnan(resid), 0, resid)
        # resid_dict
        return resid

    init_guess = [init_maxdate, 17, 1, 0, 0, 0]

    return


def make_fittedlc(fpca_f, mpfit_result, fpca_dir='', return_func=True):
    '''
    :params fpca_f: fpca filter, options: vague, B, V, R, I
    :params fpca_paramters: list of length 6
    :params covar: covariance matrix OR list of paramter errors
    :params fpca_dir: directory where fpca templates are at
    :params return_func: if true, return function, if not, return grids

    return function or grids depending on return_func
    '''

    return lambda x: x**2


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


def get_modelval(date, theta):
    '''
    Calculates the linear combination of coefficients and PC vectors at every point in the grid.
    Given a phase and parameters, returns the corresponding fitted magnitude.
    :params phase:
    :params theta:
    '''
    if not isinstance(date, float):
        date = np.array(date)
    tmax, mmax, a1, a2, a3, a4 = theta
    coeffs = [1, a1, a2, a3, a4]
    mphase = date-tmax
    y = mmax + (np.dot(coeffs,
                       np.array([fbasis(mphase) for fbasis in basis])))

    return y


def meritfunc(theta, input_data):
    '''
    Chisq to be minimized.
    :m number of samples (len(data)):
    :n number of parameters (len(theta)):
    '''
    x, y, ey = map(lambda xx: np.array(
        input_data[xx]), ['date', 'mag', 'emag'])
    y_model = get_modelval(x, theta)

    resid = (y - y_model)/ey
    resid = np.where(np.isnan(resid), 0, resid)
    # resid_dict
    return resid


theta = [0, 0, 10, 0, 0, 0]
x = np.arange(-10, 60, 1)
# print(get_modelval([0, 1], theta))
# plt.scatter(x, get_modelval(x, theta))
# plt.show()

a = pd.read_csv('02boPriv.csv')
# print(a)
data = {'date': np.array(a['MJD_OBS']) - a['MJD_OBS'].tolist()[np.argmin(np.array(a['MAG']))], 'mag': np.array(
    a['MAG']), 'emag': np.array(a['eMAG'])}
print(data['date'])
print(meritfunc(theta, data))
