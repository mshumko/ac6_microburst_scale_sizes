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
dataset_key = 'h_chorus_bwgt10'

x = np.arange(301)*50 + 25 # 50 km separation bin labels
y = np.arange(81)/40.-0.9999 # cross-correlation labels from -1 to 1.

data_array = d[dataset_key].copy()
# Estimate the probability density in each separation bin.
data_array_density = data_array/np.nansum(data_array, axis=0)
# print(data_array_density, data_array_density.shape, data_array_density[:, 1])

for i in range(20):
    plt.plot(y, data_array_density[:, i])
plt.show()

def gaus(x,a,x0,sigma):
    """ 
    Gaussian profile to fit.
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

# Now fit the profile with Gaussians.
p = np.nan*np.zeros((x.shape[0], 3)) # Array to store all Gaus parameters
for i, row in enumerate(data_array_density.T):
    #print(row, row == 'nan')
    if np.any(np.isnan(row)):
        print(i)
        continue
    p[i, :], _ = scipy.optimize.curve_fit(gaus, y, row, p0=[0.05, 0.7, 0.2])


# Calculate the median and 95% confidence interval of the density
#detection_stats = np.nanpercentile(data_array_density, [2.5, 50, 97.5], axis=0)
# detection_stats = np.nan*np.zeros((3, data_array_density.shape[1]))
# for i in range(data_array_density.shape[1]):
#     v = np.concatenate([y[i]*np.ones(int(v_i)) for i, v_i in enumerate(data_array[:, i])])
#     if len(v):
#         detection_stats[:, i] = np.percentile(v, [2.5, 50, 97.5])
#         #detection_stats[1, i] = np.median(v)
print(p)
# Colormaps attempted Greens, Plasma
pcolormesh = plt.pcolormesh(x, y, data_array_density, vmax=0.1, vmin=0.01, 
    norm=matplotlib.colors.LogNorm(), cmap=plt.get_cmap("Greens"))
plt.colorbar(pcolormesh, label='Coincidence probability')
#plt.errorbar(x, p[:, 1], yerr='none', fmt='ko')
plt.scatter(x-(x[1]-x[0])/2, p[:,1], color='k', s=3)
# plt.errorbar(x, detection_stats[1, :], errorevery=1, c='k', fmt='.')
            # yerr=[detection_stats[1, :]-detection_stats[0, :], 
            #     detection_stats[-1, :]-detection_stats[1, :]])
plt.title(dataset_key); plt.xlabel(r'THEMIS $|\Delta r|$ [km]'); plt.ylabel('correlation')
plt.xscale('log')
plt.xlim(90, 2E3)
plt.ylim(0, 1)
plt.show()