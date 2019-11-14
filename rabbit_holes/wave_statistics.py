# This script loads and plots the Agapitov et al., 2018 THEMIS chorus
# dataset using scipy.io.readsav() and looks at the distribution of
# chorus correlations as a function of THEMIS separation. Also calculates
# the median and 95% confidence interval by undoing Oleksiy's frequency 
# tables.

"""
*** DATA DESCRIPTION FROM OLEKSIY ***
the data is enclosed as arrays 301 x 81, meaning 301 intervals with 50 km 
sampling -> 0-15000 km; and 81 intervals with 0.025 step from -1 to 1 for 
the correlation coefficient.

x = findgen(301)*50.+25
y = findgen(81)/40.-0.9999

both files contain the same-named variables and I don't remember if they 
differ ) I suggest to use the vars from the corresponding file - I did so 
for my figures. You may play with the correlation coefficient threshold. 
Actually, I think we are good enough already )

Best Regards,
Oleksiy
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import scipy.io
import scipy.optimize
import os

f_dir = '/home/mike/research/ac6_microburst_scale_sizes/data'
# Filename options are: JGR2018_F_Fig4_Hchorus_BW_gt_10.sav or 
# JGR2018_F_Fig4_Hchorus_BW_gt_lt_10.sav
f_name = 'JGR2018_F_Fig4_Hchorus_BW_gt_lt_10.sav' 
d = scipy.io.readsav(os.path.join(f_dir, f_name))
dataset_key = 'h_chorus_bwlt10'

x = np.arange(301)*50 + 25 # 50 km separation bin labels
y = np.arange(81)/40.-0.9999 # cross-correlation labels from -1 to 1.

MAX_SEP = 2000
max_idx = np.where(x > MAX_SEP)[0][0]

data_array = d[dataset_key].copy()
# Estimate the probability density in each separation bin.
data_array_density = data_array/np.nansum(data_array, axis=0)

# for i in range(max_idx):
#     plt.plot(y, data_array_density[:, i])
# plt.show()

def gaus(x,a,x0,sigma):
    """ 
    Gaussian profile to fit.
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

# Now fit the profile with Gaussians.
p = np.nan*np.zeros((max_idx, 3)) # Array to store all Gaus parameters
for i, row in enumerate(data_array_density[:, :max_idx].T):
    # Don't try to fit NaNs
    if np.any(np.isnan(row)):
        continue
    p[i, :], _ = scipy.optimize.curve_fit(gaus, y, row, p0=[0.05, 0.7, 0.2])

print(p)


######################## PLOTTING ##########################
# Colormaps attempted Greens, Plasma
pcolormesh = plt.pcolormesh(x, y, data_array_density, vmax=0.1, vmin=0.01, 
    norm=matplotlib.colors.LogNorm(), cmap=plt.get_cmap("Greens"))
plt.colorbar(pcolormesh, label='Chorus coincidence probability')
#plt.errorbar(x, p[:, 1], yerr='none', fmt='ko')
plt.scatter(x[:max_idx]-(x[1]-x[0])/2, p[:,1], color='k', s=3)

if 'bwgt' in dataset_key:
    plt.title(r'$B_w > $ 10 pT'); 
if 'bwlt' in dataset_key:
    plt.title(r'$B_w < $ 10 pT'); 

plt.xlabel(r'THEMIS $|\Delta r|$ [km]'); plt.ylabel('correlation')
# plt.xscale('log')
plt.xlim(100, MAX_SEP)
plt.ylim(0, 1)
plt.show()