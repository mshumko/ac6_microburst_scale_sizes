# This script plots the distribution of AE vs microburst detections.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

catalog_path = (f'/home/mike/research/'
                    'ac6_microburst_scale_sizes'
                    '/data/coincident_microbursts_catalogues/'
                    'AC6_coincident_microbursts_sorted'
                    f'_v{6}.txt')
cat = pd.read_csv(catalog_path)

# plt.scatter(cat['Dist_Total'], cat['AE'])
# plt.show()

ae_bins_nt = np.arange(0, 1000, 50)
s_bins_km = np.arange(0, 101, 10)
H = np.nan*np.zeros((len(s_bins_km), len(ae_bins_nt)-1))

for i, (s_bin_i, s_bin_f) in enumerate(zip(s_bins_km[:-1], s_bins_km[1:])):
    df_flt = cat[(cat['Dist_Total'] > s_bin_i) & (cat['Dist_Total'] <= s_bin_f)]
    H[i, :], _ = np.histogram(df_flt.AE, bins=ae_bins_nt, density=True)
    # plt.step(ae_bins_nt[:-1], H[i, :], where='post', label=f'{s_bin_i}-{s_bin_f} [km]')
    plt.plot(ae_bins_nt[:-1], H[i, :], label=f'{s_bin_i}-{s_bin_f} [km]')

plt.legend()
plt.show()