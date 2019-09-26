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
import os

f_dir = '/home/mike/research/ac6_microburst_scale_sizes/data'
# Filename options are: JGR2018_F_Fig4_Hchorus_BW_gt_10.sav or 
# JGR2018_F_Fig4_Hchorus_BW_gt_lt_10.sav
f_name = 'JGR2018_F_Fig4_Hchorus_BW_gt_lt_10.sav' 
d = scipy.io.readsav(os.path.join(f_dir, f_name))
dataset_key = 'h_chorus_bwlt10'

x = np.arange(301)*50 + 25 # 50 km separation bin labels
y = np.arange(81)/40.-0.9999 # cross-correlation labels from -1 to 1.

# Convert bins with 0 detections to np.nan
data_array = d[dataset_key].copy()

# Normlize the h_chorus_bwgt10 histrogram in the separation dimension
max_det = np.nansum(d[dataset_key], axis=0) # max detections in each separation bin.
np.tile(max_det, 81).reshape(*data_array.shape)
data_array_density = data_array/np.tile(max_det, 81).reshape(*data_array.shape)

# Calculate the median and 95% confidence interval of the density
#detection_stats = np.nanpercentile(data_array_density, [2.5, 50, 97.5], axis=0)
detection_stats = np.nan*np.zeros((3, data_array_density.shape[1]))
for i in range(data_array_density.shape[1]):
    v = np.concatenate([y[i]*np.ones(int(v_i)) for i, v_i in enumerate(data_array[:, i])])
    if len(v):
        detection_stats[:, i] = np.percentile(v, [2.5, 50, 97.5])
        #detection_stats[1, i] = np.median(v)

# Colormaps attempted Greens, Plasma
p = plt.pcolormesh(x, y, data_array_density, vmax=0.1, vmin=0.01, 
    norm=matplotlib.colors.LogNorm(), cmap=plt.get_cmap("Greens"))
plt.colorbar(p, label='Coincidence probability')
plt.errorbar(x, detection_stats[1, :], errorevery=1, c='k')
            # yerr=[detection_stats[1, :]-detection_stats[0, :], 
            #     detection_stats[-1, :]-detection_stats[1, :]])
plt.title(dataset_key); plt.xlabel(r'THEMIS $|\Delta r|$ [km]'); plt.ylabel('correlation')
plt.xscale('log')
plt.xlim(100, 1E4)
plt.ylim(0, 1)
plt.show()