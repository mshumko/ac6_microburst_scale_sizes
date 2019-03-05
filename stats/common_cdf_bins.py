# This script finds the 10 most common bins for the microburst CDF as a function of separation.
import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser
import itertools
import os

import pandas as pd

# Load microburst catalog
version = 5
catalog_path = ('/home/mike/research/ac6_microburst_scale_sizes/'
                'data/microburst_catalogues'
                '/AC6A_microbursts_v{}.txt'.format(version))
# converters = {
#             0:dateutil.parser.parse, 
#             -1:dateutil.parser.parse, 
#             -2:dateutil.parser.parse
#             }
data = pd.read_csv(catalog_path)

# Binned counts directory
bin_counts_dir = ('/home/mike/research/ac6_microburst_scale_sizes/'
            'data/binned_counts')

### FILTERS ###
# Filter out detections outside of the outer radiation belt
# and above the US.
# L filter
data = data[(data['Lm_OPQ'] > 4) & (data['Lm_OPQ'] < 8)]
# USA filter
data = data[
    ((data['lon'] > -60) | (data['lon'] < -140)) |
    ((data['lat'] > 70) | (data['lat'] < 15))
    ]
# SAA filter
#data = data[
#    ((data['lon'] > 30)  | (data['lon'] < -116)) |
#    ((data['lat'] < -90) | (data['lat'] > 0))
#    ]
# Filter by the number of standrad deviations a peak is 
# above a 10% baseline.
#data = data[data['peak_std'] > 2]
# Filter out ambigious CCs with curtains
# data = data[data['time_cc'] > data['space_cc']+curtain_thresh]

### Create all of the L-MLT-AE bins to iterate over. ###
dL = 1
dMLT = 1
dAE = 100
L_bins=np.arange(4, 8, dL) 
MLT_bins=np.arange(0, 24, dMLT)
AE_bins=np.arange(0, 600, dAE)

# Create all possible combinations of L, MLT and AE
LL, MLTMLT, AEAE = np.meshgrid(L_bins, MLT_bins, AE_bins)
bin_vals_tuple = zip(LL.flatten(), MLTMLT.flatten(), AEAE.flatten())
len_bins = len(LL.flatten())

# Define arrays to save.
N_bursts = np.nan*np.zeros(len_bins)
N_bin_samples = np.nan*np.zeros(len_bins)
N_burst_rate = np.nan*np.zeros(len_bins)

def index_to_name(L, MLT, AE):
    s = '{}_L_{}_{}_MLT_{}_{}_AE_{}'.format(
        L, L+dL, MLT, MLT+dMLT, AE, AE+dAE)
    return s

### BIN LOOP ###
for i, (Li, MLTi, AEi) in enumerate(bin_vals_tuple):
    filtered_catalog = data[
        (data['Lm_OPQ'] > Li) & (data['Lm_OPQ'] < Li+dL) &
        (data['MLT_OPQ'] > MLTi) & (data['MLT_OPQ'] < MLTi+dMLT) &
        (data['AE'] > AEi) & (data['AE'] < AEi+dAE)
        ]
    # Number of microbursts detected in that L-MLT-AE bin.
    N_bursts[i] = filtered_catalog.shape[0] 

    # If enough microbursts were detected in that bin.
    if N_bursts[i] >= 100: 
        # Load the current count bin
        s = index_to_name(Li, MLTi, AEi)
        binName = 'AC6_counts_' + s + '.csv'
        binPath = os.path.join(bin_counts_dir, binName)
        binDf = pd.read_csv(binPath)
        N_bin_samples[i] = binDf.shape[0]
        N_burst_rate[i] = round(100*N_bursts[i]/N_bin_samples[i], 2)

# Put all of the data into one array.
d = np.stack(
    (LL.flatten(), MLTMLT.flatten(), AEAE.flatten(), 
     N_bursts, N_bin_samples, N_burst_rate),
     axis=1)
# Load data into a DataFrame.
columns = ['L', 'MLT', 'AE', 'bursts', 'samples', 'rate']
df = pd.DataFrame(data=d, dtype=float, columns=columns)
# Remove nan values
df = df.dropna()
# Sort by the microburst rate.
df = df.sort_values(by='rate', ascending=False)
# write to csv file.
df.to_csv('most_common_l_mlt_ae_bins.csv', index=False)