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

snlcpy accepts several input types for data. For date, magnitude, and magnitude error, you may use three lists, three numpy arrays, an astropy table, a pandas dataframe, or a dictionary. Magnitude error is an optional input. If you include magnitude error, the FPCA fit will be weighted by these errors. 

snlcpy only fits one photometric bandpass at a time. Available templates are Johnson *B, V, R, I* filters and one 'vague' template, which may be used for any bandpass. If a Johnson filter is not specified in fit_params(), the program will default to 'vague'. 

### Example
#### Using a pandas dataframe
```
sn = pd.read_csv('SN2011fe.csv')
B_data = sn[sn['filt'] == 'B']
INCOMPLETE EXAMPLE

```
#### Using an astropy table
```
INCOMPLETE EXAMPLE
```
#### Using three lists or numpy arrays
```
sn = pd.read_csv('SN2011fe.csv')
B_data = sn[sn['filt'] == 'B']
date = np.array(B_data['JD'])
mag = np.array(B_data['mag'])
emag = np.array(B_data['mag_err'])
res = snlcpy.fit_pcparams([date,mag,emag],'B')
snlcpy.plot_lc([date,mag,emag],res,fpca_f='B',input_func=True)
```
#### Using a dictionary
```
sn = pd.read_csv('SN2011fe.csv')
B_data = sn[sn['filt'] == 'B']
B_data_sub = B_data[['JD', 'mag', 'mag_err']]
B_data_dict = B_data_sub.to_dict('list')
res = snlcpy.fit_pcparams(B_data_dict,'B')
snlcpy.plot_lc(B_data_dict,res,fpca_f='B',input_func=True)
```