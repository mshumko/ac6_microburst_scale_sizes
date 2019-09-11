# This script reads in Oleksiy's data and makes a plot of 
# the probability that a highly correlated chorus wave 
# (correlation > 0.8) is observed at each separation.

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors
import scipy.io
import os

f_dir = '/home/mike/research/ac6_microburst_scale_sizes/data'
# Filename options are: JGR2018_F_Fig4_Hchorus_BW_gt_10.sav or 
# JGR2018_F_Fig4_Hchorus_BW_gt_lt_10.sav
f_name = 'JGR2018_F_Fig4_Hchorus_BW_gt_lt_10.sav' 
wave_data = scipy.io.readsav(os.path.join(f_dir, f_name))

# Specify a few helper variables.
gt_key = 'h_chorus_bwgt10'
lt_key = 'h_chorus_bwlt10'
n_c, n_s = wave_data[gt_key].shape
cc_thresh = 0.8

# Oleksiy's bins
x = np.arange(n_s)*50 + 25 # 50 km separation bin labels
y = np.arange(n_c)/40.-0.9999 # cross-correlation labels from -1 to 1.

idc = np.where(y >= cc_thresh)[0][0] # Index of cc_thresh correlation bin edge

gt_fraction = np.sum(wave_data[gt_key][idc:, :], axis=0)/np.sum(wave_data[gt_key], axis=0)
lt_fraction = np.sum(wave_data[lt_key][idc:, :], axis=0)/np.sum(wave_data[lt_key], axis=0)

# Normalize to max of 1
# gt_fraction /= np.nanmax(gt_fraction)
# lt_fraction /= np.nanmax(lt_fraction)

gt_fraction_err = gt_fraction*(1-gt_fraction)/np.sqrt(np.sum(wave_data[gt_key][idc:, :], axis=0))
lt_fraction_err = lt_fraction*(1-lt_fraction)/np.sqrt(np.sum(wave_data[lt_key][idc:, :], axis=0))


plt.errorbar(x, gt_fraction, yerr=gt_fraction_err, c='r', label='bw > 10 pT')
plt.errorbar(x, lt_fraction, yerr=lt_fraction_err, c='b', label='bw < 10 pT')
plt.legend()
plt.xlim(100, 1E3)
plt.show()