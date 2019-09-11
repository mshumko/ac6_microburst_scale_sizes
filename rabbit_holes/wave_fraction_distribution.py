# This script reads in Oleksiy's data and makes a plot of 
# the probability that a highly correlated chorus wave 
# (correlation > 0.8) is observed at each separation.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.io
import os

f_dir = '/home/mike/research/ac6_microburst_scale_sizes/data'
# Filename options are: JGR2018_F_Fig4_Hchorus_BW_gt_10.sav or 
# JGR2018_F_Fig4_Hchorus_BW_gt_lt_10.sav
f_name = 'JGR2018_F_Fig4_Hchorus_BW_gt_lt_10.sav' 
wave_data = scipy.io.readsav(os.path.join(f_dir, f_name))

# Read the equatorial microburst fraction
file_dir = '/home/mike/research/ac6_microburst_scale_sizes/data'
file_name = 'equatorial_microburst_fraction.csv'
file_path = os.path.join(file_dir, file_name)
microburst_fraction = pd.read_csv(file_path, index_col='s')

# Specify a few helper variables.
gt_key = 'h_chorus_bwgt10'
lt_key = 'h_chorus_bwlt10'
n_c, n_s = wave_data[gt_key].shape
cc_thresh = 0.7

# Oleksiy's bins
x = np.arange(n_s)*50 + 25 # 50 km separation bin labels
y = np.arange(n_c)/40.-0.9999 # cross-correlation labels from -1 to 1.

idc = np.where(y >= cc_thresh)[0][0] # Index of cc_thresh correlation bin edge

gt_fraction = np.sum(wave_data[gt_key][idc:, :], axis=0)/np.sum(wave_data[gt_key], axis=0)
lt_fraction = np.sum(wave_data[lt_key][idc:, :], axis=0)/np.sum(wave_data[lt_key], axis=0)

# Normalize to max of 1
gt_fraction /= np.nanmax(gt_fraction)
lt_fraction /= np.nanmax(lt_fraction)
microburst_fraction.f /= np.max(microburst_fraction.f)
microburst_fraction.f_err /= np.max(microburst_fraction.f)

# Calculate error
gt_fraction_err = gt_fraction*(1-gt_fraction)/np.sqrt(np.sum(wave_data[gt_key][idc:, :], axis=0))
lt_fraction_err = lt_fraction*(1-lt_fraction)/np.sqrt(np.sum(wave_data[lt_key][idc:, :], axis=0))


plt.errorbar(x, gt_fraction, yerr=gt_fraction_err, c='r', label='bw > 10 pT')
plt.errorbar(x, lt_fraction, yerr=lt_fraction_err, c='b', label='bw < 10 pT')
plt.errorbar(microburst_fraction.index, microburst_fraction.f, 
            yerr=microburst_fraction.f_err, c='k', label='microbursts')
plt.legend()
plt.xlim(100, 1E3)
plt.ylim(0, 1)
plt.title(f'Coincident probability of chorus waves vs microbursts\ncc_thresh = {cc_thresh}')
plt.xlabel('Equatorial separation [km]')
plt.ylabel('Normalized coincident probability')
plt.show()