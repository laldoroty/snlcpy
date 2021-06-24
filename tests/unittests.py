import pytest
from astropy.table import Table
from snlcpy import fpca
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
abspath = os.path.dirname(__file__)

filters = ['B','V','R','I','vague']
test_grid = np.arange(-10,40.1,0.1)
testsn = pd.read_csv(os.path.join(abspath, '02boPriv.csv'))

def test_getpctemplates():
    '''
    Tests that the FPCA templates are pulled.
    '''
    for filt in filters:
        PCs = fpca.get_pctemplates(filt)
        assert len(PCs) == 5

def test_make_fittedlc():
    for filt in filters:
        fpca.make_fittedlc(filt,)

# def test_get_modelval():
#     theta = np.array([0,0,0,0,0,0])
#     for filt in filters:
#         PCs = fpca.get_pctemplates(filt)
#         testlc = PCs[0](test_grid)

#         assert np.allclose(testlc,fpca.get_modelval(test_grid,theta))
#         print(filt, 'passes')




# def test_fit_pcparams():
#     for filt in filters:
#         PCs = fpca.get_pctemplates(filt)
#         test_grid = np.arange(-10,40,0.1)
#         fpca.fit_pcparams(test_grid,PCs[0],fpca_f=filt,
#             init_guess=[0,0,0,0,0,0],components=4)

# # test_getpctemplates()
# test_fit_pcparams()
# test_get_modelval()