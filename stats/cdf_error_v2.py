# PDF and CDF counting error calculation using Monte Carlo methods
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

NORM = False
s_bins_km = np.arange(0, 100, 5)

# Load microburst catalog
catalog_version = 6
catalog_path = (f'/home/mike/research/'
                    'ac6_microburst_scale_sizes'
                    '/data/coincident_microbursts_catalogues/'
                    'AC6_coincident_microbursts_sorted'
                    f'_v{catalog_version}.txt')
cat = pd.read_csv(catalog_path)

if NORM:
    # Load normalization dataset
    norm_path = ('/home/mike/research/ac6_microburst'
                '_scale_sizes/data/norm/ac6_norm_all_cdf.csv')
    norm_data = pd.read_csv(norm_path, index_col=0)
    min_s = s_bins_km[0]
    max_s = s_bins_km[-2]
    norm = norm_data.Seconds.max()/norm_data
    norm = norm.loc[min_s:max_s, 'Seconds']

# Histogram the catalog data
H, _ = np.histogram(cat['Dist_Total'], bins=s_bins_km)

if NORM:
    H = H * norm

# Now wiggle the number of detections in each bin assuming a 
# Poisson distribution.
N = 10000
pdf_trials = np.nan*np.ones((len(H), N), dtype=int)
for i, H_i in enumerate(H):
    pdf_trials[i, :] = np.random.poisson(H_i, size=N)

# Now calculate the CDF for all trials.
cdf  = 100*np.array([np.sum(H[j:]) for j in range(len(H))])/np.sum(H)
cdf_trials = np.nan*np.ones_like(pdf_trials)
for i, pdf_i in enumerate(pdf_trials.T):
    cdf_trials[:, i] = 100*np.array([np.sum(pdf_i[j:]) for j in range(len(pdf_i))])/np.sum(pdf_i)

# Now calculate the Q1-3 range for the PDF.
pdf_q = np.percentile(pdf_trials, [25, 75], axis=1)
cdf_q = np.percentile(cdf_trials, [25, 75], axis=1)

# Plot everything
fig, ax = plt.subplots(2, sharex=True)
ax[0].fill_between(s_bins_km[:-1], pdf_q[0, :], pdf_q[1, :], step='post', 
                    label='Monte Carlo Q1-Q3')
ax[0].step(s_bins_km[:-1], H, c='r', where='post', label='Cataloged microbursts')

ax[1].step(s_bins_km[:-1], cdf, c='r', where='post')
ax[1].fill_between(s_bins_km[:-1], cdf_q[0, :], cdf_q[1, :], step='post', 
                    label='Monte Carlo Q1-Q3')

ax[0].set_title(f'AC6 Microburst Distribution\nPoisson Error Monte Carlo | Normalized={NORM}')
ax[1].set_xlabel('AC6 separation [km]')
ax[0].set_ylabel('Number of microbursts')
ax[1].set_xlim(left=0)
ax[1].set_ylabel('F(s) [%]')
ax[0].legend()
plt.show()