import pytest
from astropy.table import Table
from snlcpy import fpca

def test_getpctemplates():
    '''
    Tests that the correct FPCA templates are pulled.
    '''
    filters = ['B','V','R','I']
    print(filters)

    for filt in filters:
        PCs = fpca.get_pctemplates(filt)
        print(filt)

test_getpctemplates()