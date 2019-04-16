# This script calls the 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import csv

import ac6_microburst_scale_sizes.microburst_detection.replace_error_sep_lags as replace_error_sep_lags

r = replace_error_sep_lags.ReplaceErrorVals('/home/mike/research/ac6/AC6_Separation.csv', None)
r.loadSeprationFile()

# Now load the L-MLT normalization files.
with open('/home/mike/research/ac6_microburst_scale_sizes/data/norm/ac6_L_MLT_bins.csv') as f:
    keys = next(f).rstrip().split(',')
    bins = {}
    for key in keys:
        bins[key] = next(f).rstrip().split(',')
        bins[key] = list(map(float, bins[key]))
with open('/home/mike/research/ac6_microburst_scale_sizes/data/norm/ac6_L_MLT_norm.csv') as f:
    reader = csv.reader(f)
    next(reader) # skip header
    norm = 10*np.array(list(reader)).astype(float) # Convert to number of samples.

# Plot the lifetime AC-6 separation
fig, ax = plt.subplots(2)
ax[0].plot(r.sepDict['Date/Time'], np.abs(r.sepDict['In-Track Separation [km]']))
ax[0].set_xlabel('Date [YYYY-MM]')
ax[0].set_ylabel('Separation [km]')
ax[0].set_title('AC-6 In-Track Separation')
ax[0].set_yscale('log')
ax[0].set_ylim(bottom=1)
ax[0].set_xticks(pd.date_range(start='6/1/2014', end='9/1/2017', freq='2Q'))
ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))

# L shell filter for the L-MLT plot
L_lower = 4
idL = np.where(np.array(bins['Lm_OPQ']) >= L_lower)[0][0]
p = ax[1].pcolormesh(bins['MLT_OPQ'], bins['Lm_OPQ'][idL:], norm[idL:, :]/1E5, cmap='Greens')
plt.colorbar(p, ax=ax[1], label=r'Samples x $10^5$')
ax[1].set_xlabel('MLT')
ax[1].set_ylabel('L')

plt.tight_layout()
plt.show()