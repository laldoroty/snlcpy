import pytest
from astropy.table import Table
from snlcpy import fpca

def test_getpctemplates():
    '''
    Tests that the FPCA templates are pulled.
    '''
    filters = ['B','V','R','I','vague']
    for filt in filters:
        PCs = fpca.get_pctemplates(filt)
        assert len(PCs) == 5

test_getpctemplates()