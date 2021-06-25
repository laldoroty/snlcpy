# snlcpy
Python package that uses the results from [He et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...857..110H/abstract) in order to fit light curves to SNe Ia. Please cite this paper, as well as [ZENODO LINK], when using this code. 

## Installation
1. **Before installing snlcpy**, you need [pycmpfit](https://github.com/cosmonaut/pycmpfit). pycmpfit has a dependency on cython. Use the following code to download and install pycmpfit:

```
git clone https://github.com/cosmonaut/pycmpfit
python setup.py build_ext
```

## Using snlcpy

snlcpy accepts several input types for data. For date, magnitude, and magnitude error, you may use three lists, three numpy arrays, an astropy table, or a dictionary. Magnitude error is an optional input. If you include magnitude error, the FPCA fit will be weighted by these errors. 

snlcpy only fits one photometric bandpass at a time. Available templates are Johnson *B, V, R, I* filters and one 'vague' template, which may be used for any bandpass. If a Johnson filter is not specified in fit_params(), the program wil default to 'vague'. 

### Example

```
sn = pd.read_csv('SN2011fe.csv')
data = {'date': epoch, 'mag': mag, 'emag': emag}
res = fit_params(data,'B',components=2)
lc = make_fittedlc('B', res)
```
