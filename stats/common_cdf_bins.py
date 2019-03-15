# This script finds the 10 most common bins for the microburst CDF as a function of separation.
import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser
import itertools
import os

import pandas as pd

# Load microburst catalog
version = 5
catalog_dir = ('/home/mike/research/ac6_microburst_scale_sizes/'
                'data/microburst_catalogues')

catalog_name = f'/AC6A_microbursts_v{version}.txt'
# catalog_dir = ('/home/mike/research/ac6_microburst_scale_sizes/'
#                             'data/coincident_microbursts_catalogues')
# catalog_name = f'AC6_coincident_microbursts_sorted_v{version}.txt'
data = pd.read_csv(os.path.join(catalog_dir, catalog_name))

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

# Number of highest rate bins to save.
N_SAVE_BINS = 5 

### Create all of the L-MLT-AE bins to iterate over. ###
dL = 1
dMLT = 1
dAE = 100
d_sep = 5

L_bins=np.arange(4, 8, dL) 
MLT_bins=np.arange(0, 24, dMLT)
AE_bins=np.arange(0, 600, dAE)
D_bins = np.arange(0, 200, d_sep)

# Create all possible combinations of L, MLT and AE
LL, MLTMLT, AEAE = np.meshgrid(L_bins, MLT_bins, AE_bins)
bin_vals_tuple = list(zip(LL.flatten(), MLTMLT.flatten(), AEAE.flatten()))
len_count_bins = len(LL.flatten())

# Define a dictionary of len_count_bins x 5 arrays for each separation
# bin. First three columns are L, MLT, and AE lower edge bins, and then
# the number of detections, and lastly the microburst rate.
bin_burst_num = {d:np.nan*np.ones((len_count_bins, 5), dtype='float')
                    for d in D_bins}
count_bin_samples = np.nan*np.ones(len_count_bins, dtype=int)

def index_to_name(L, MLT, AE):
    """ 
    Helper function to convert from L, MLT, and AE
    values to a string to load count bin file.
    """
    s = '{}_L_{}_{}_MLT_{}_{}_AE_{}'.format(
        L, L+dL, MLT, MLT+dMLT, AE, AE+dAE)
    return s

### COUNT BIN SAMPLE SIZE LOOP ###
for i, (Li, MLTi, AEi) in enumerate(bin_vals_tuple):
    s = index_to_name(Li, MLTi, AEi)
    binName = 'AC6_counts_' + s + '.csv'
    binPath = os.path.join(bin_counts_dir, binName)
    # If that bin exists (some bins are rare or non-existant!)
    if os.path.exists(binPath):
        binDf = pd.read_csv(binPath)
        # Save the number of data points in that bin.
        count_bin_samples[i] = binDf.shape[0]

### BIN LOOP ###
# Loop over the separation bins
for d_i, d in enumerate(D_bins): 
    # Loop over the count bins.
    for i, (Li, MLTi, AEi) in enumerate(bin_vals_tuple):
        filtered_catalog = data[
            # L shell filter
            (data['Lm_OPQ'] > Li) & (data['Lm_OPQ'] < Li+dL) &
            # MLT filter
            (data['MLT_OPQ'] > MLTi) & (data['MLT_OPQ'] < MLTi+dMLT) &
            # AE filter
            (data['AE'] > AEi) & (data['AE'] < AEi+dAE) & 
            # Total distance filter
            (data['Dist_Total'] > d)
            ]
        # Number of microbursts detected in that L-MLT-AE bin.
        bin_burst_num[d][i, 0:3] = [Li, MLTi, AEi]
        bin_burst_num[d][i, 3] = filtered_catalog.shape[0] 

        # If enough microbursts were detected in that bin.
        if filtered_catalog.shape[0] >= 100: 
            # Load the current count bin
            s = index_to_name(Li, MLTi, AEi)
            binName = 'AC6_counts_' + s + '.csv'
            binPath = os.path.join(bin_counts_dir, binName)
            binDf = pd.read_csv(binPath)
            # Save the burst rate in units of microbursts/seconds
            bin_burst_num[d][i, 4] = round(10*bin_burst_num[d][i, 0]/count_bin_samples[i], 4)

# Organize this data into a dictionary of DataFrames
columns = ['L', 'MLT', 'AE', 'burst_num', 'burst_rate']
df_dict = {key:pd.DataFrame(data=bin_burst_num[key], columns=columns) 
                    for key in bin_burst_num.keys()}
# Loop over each cumulative separation bin and sort by the 
# microburst rates. 
for key, df in df_dict.items():
    # Drop nan rows
    df_dict[key] = df_dict[key].dropna()
    # Sort by microburst rates
    df_dict[key] = df_dict[key].sort_values('burst_rate', ascending=False)
    # Keep only N_SAVE_BINS count bins with highest occurance rates
    df_dict[key] = df_dict[key].iloc[:N_SAVE_BINS]
    # Insert a separation row.
    df_dict[key].insert(0, 'd', key)

    # If first time around write the header. Otherwise do not.
    if os.path.exists('most_common_counts_bins.csv'):
        header=False
    else:
        header=True
    # Save to a csv file.
    with open('most_common_counts_bins.csv', 'a') as f:
        df_dict[key].to_csv(f, header=header, index=False)

# # Put all of the data into one array.
# d = np.stack(
#     (LL.flatten(), MLTMLT.flatten(), AEAE.flatten(), 
#      N_bursts, N_bin_samples, N_burst_rate),
#      axis=1)
# # Load data into a DataFrame.
# columns = ['L', 'MLT', 'AE', 'bursts', 'samples', 'rate']
# df = pd.DataFrame(data=d, dtype=float, columns=columns)
# # Remove nan values
# df = df.dropna()
# # Sort by the microburst rate.
# df = df.sort_values(by='rate', ascending=False)
# # write to csv file.
# df.to_csv('most_common_l_mlt_ae_bins.csv', index=False)