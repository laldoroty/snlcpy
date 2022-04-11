# snlcpy
Python package that uses the results from [He et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...857..110H/abstract) in order to fit light curves to SNe Ia. Please cite this paper, as well as [ZENODO LINK], when using this code. 

Package by [Lauren Aldoroty](https://laldoroty.github.io), Jiawen Yang 2021.

Package currently incomplete. 

## Installation

```
git clone https://github.com/laldoroty/snlcpy
```

Then, `cd` into the `snlcpy` directory.

```
python setup.py install
```

## Using snlcpy

snlcpy accepts three lists or three numpy arrays as data input. Magnitude error is an optional input. If you include magnitude error, the FPCA fit will be weighted by these errors. 

snlcpy only fits one photometric bandpass at a time. Available templates are Johnson *B, V, R, I* filters and one 'vague' template, which may be used for any bandpass. If a Johnson filter is not specified in fit_params(), the program will default to 'vague'. 

### Example
```
sn = pd.read_csv('SN2011fe.csv')
B_data = sn[sn['filt'] == 'B']
date = np.array(B_data['JD'])
mag = np.array(B_data['mag'])
emag = np.array(B_data['mag_err'])
res = snlcpy.fit_pcparams(date,mag,emag,'B')
snlcpy.plot_lc(date,mag,res,fpca_f='B',emag=emag,input_func=True)
```

The fit result (i.e., PC coefficients) from this example can be accessed with:

```
print(res['params']))
```
