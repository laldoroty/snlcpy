
import numpy as np


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
    ## Read in data
    if len(data) == 3:
        date, mag, emag = [np.array(ii) for ii in data]
    if len(data) == 2:
        date, mag = [np.array(ii) for ii in data]
        # if no magnitude error is provided, then use constant 1 for all epochs
        emag = np.array([1]*len(date))

    if len(data) == 1 and isinstance(data[0], dict):
        data = data[0]
        date, mag = map(lambda x: np.array(data[x]), ['date', 'mag'])
        if 'emag' not in data:
            emag = np.array([1]*len(date))
        else:
            emag = np.array(data['emag'])

    # Initial guess
    init_maxdate = date[np.argmin(mag)]

    print(date, mag, emag)
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


fit_pcparams({'date': [1, 2, 3], 'mag': [3, 4, 5], 'emag': [4, 6, 3]})
