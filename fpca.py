import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d

# Use for testing:
from numpy.core.fromnumeric import shape


def fit_pcparams(*data, mag=None, magerr=None, fpca_f='vague',
                 init_guess=None, components=2, fpca_dir='',
                 penalty=True, penalty_increase=None, penalty_decrease=None):
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

    print(date, mag, emag)

    # Initial guess
    init_maxdate = date[np.argmin(mag)]
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


# test_dict = {'date': [1, 2, 3], 'mag': [3, 4, 5], 'emag': [4, 6, 3]}
# fit_pcparams(test_dict)
test_tab = Table([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                 names=('date', 'mag', 'emag'))
print('success')
fit_pcparams(test_tab)
print('success #2')