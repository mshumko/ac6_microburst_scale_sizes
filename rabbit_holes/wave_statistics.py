# This script loads and plots the Agapitov et al., 2018 THEMIS chorus
# dataset using scipy.io.readsav()

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
dataset_key = 'h_chorus_bwgt10'

x = np.arange(301)*50 + 25 # 50 km separation bin labels
y = np.arange(81)/40.-0.9999 # cross-correlation labels from -1 to 1.

# Normlize the h_chorus_bwgt10 histrogram in the separation dimension
max_det = np.sum(d[dataset_key], axis=0) # max detections in each separation bin.
np.tile(max_det, 81).reshape(*d[dataset_key].shape)
d[f'{dataset_key}_density'] = d[dataset_key]/np.tile(max_det, 81).reshape(*d[dataset_key].shape)

p = plt.pcolormesh(x, y, d[f'{dataset_key}_density'], vmax=0.1, vmin=0.01, norm=matplotlib.colors.LogNorm())
    #cmap=matplotlib.colors.Colormap('Greens'))
plt.colorbar(p)
plt.xscale('log')
plt.xlim(100, 1E4)
plt.ylim(0, 1)
plt.show()