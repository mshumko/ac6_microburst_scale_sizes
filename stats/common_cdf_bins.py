# This script finds the 10 most common bins for the microburst CDF as a function of separation.
import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser
import itertools
import os

import pandas as pd

# Load coincident microburst catalog
version = 4
catPath = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
        'coincident_microbursts_catalogues/'
       'AC6_coincident_microbursts_v{}.txt'.format(version))
converters = {
            0:dateutil.parser.parse, 
            -1:dateutil.parser.parse, 
            -2:dateutil.parser.parse
            }
data = pd.read_csv(catPath, converters=converters)

### FILTERS ###
# Filter out detections outside of the outer radiation belt
# and above the US and SAA. Also filter out detections
# which had a higher spatial_CC by a curtain_thresh ammount. 
# Lastly, filter out detections that are close to the noise
# floor (peak_std)
curtain_thresh = 0.3
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
data = data[data['time_cc'] > data['space_cc']+curtain_thresh]

### BIN THE DETECTIONS ###
# Now bin the data by L-MLT-AE and separation.
L_bins=np.arange(4, 8, 1) 
MLT_bins=np.arange(0, 15, 1)
AE_bins=np.arange(0, 600, 100)
sep_bins = np.arange(0, 100, 5)

def index_to_name(LL, MLTMLT, AEAE, i, j, k):
    s = '{}_L_{}_{}_MLT_{}_{}_AE_{}'.format(
        LL[i, j, k], LL[i, j+1, k], MLTMLT[i, j, k], MLTMLT[i+1, j, k],
        AEAE[i, j, k], AEAE[i, j, k+1])
    return s

LL, MLTMLT, AEAE = np.meshgrid(L_bins, MLT_bins, AE_bins)
N = np.nan*np.zeros_like(LL)
lL, lMLT, lAE = LL.shape

binDir = ('/home/mike/research/ac6_microburst_scale_sizes/'
            'data/binned_counts')

N_most_common = 5
# most_common_bins = np.array([[None]*(len(sep_bins)-1)]*N_most_common)
most_common_bins = np.ones((len(sep_bins), 1+N_most_common*2), dtype=object)

# Loop over the L, MLT, AE bins
for s_i, d in enumerate(sep_bins):
    N = np.nan*np.zeros(((lL-1)*(lMLT-1)*(lAE-1), 4), dtype=float)
    n = 0

    for i, j, k in itertools.product(range(lL-1), range(lMLT-1), range(lAE-1)):
        data_filtered = data[
            (data['Lm_OPQ'] > LL[i, j, k]) & (data['Lm_OPQ'] < LL[i, j+1, k]) &
            (data['MLT_OPQ'] > MLTMLT[i, j, k]) & (data['MLT_OPQ'] < MLTMLT[i+1, j, k]) &
            (data['AE'] > AEAE[i, j, k]) & (data['AE'] < AEAE[i, j, k+1]) & 
            (data['Dist_Total'] > d)
            ]
        N_microburst = data_filtered.shape[0]
        if N_microburst:
            # Load the current count bin to get microburst rate
            s = index_to_name(LL, MLTMLT, AEAE, i, j, k)
            binName = 'AC6_counts_' + s + '.csv'
            binPath = os.path.join(binDir, binName)
            binDf = pd.read_csv(binPath)
            print('N_microburst/binDf.shape[0] =', 10*100*N_microburst/binDf.shape[0])
            N[n] = [i, j, k, 10*100*N_microburst/binDf.shape[0]]
            print(N[n])
        n += 1

    df = pd.DataFrame(data=N, dtype=int)
    #print(df.nlargest(N_most_common, 3))
    nlargest = df.nlargest(N_most_common, 3)
    most_common_bins[s_i, 0] = d
    for row in range(nlargest.shape[0]):
        ii, jj, kk, nn = nlargest.iloc[row]
        most_common_bins[s_i, 1+row] = index_to_name(LL, MLTMLT, AEAE, ii, jj, kk)
        most_common_bins[s_i, 1+N_most_common+row] = nn

# Save the data to file.
df_most_common = pd.DataFrame(data=most_common_bins)
header = ['lower_total_dist'] + ['{}_common'.format(i+1) for i in range(N_most_common)] + ['{}_rate_%'.format(i+1) for i in range(N_most_common)]
df_most_common.to_csv('most_common_l_mlt_ae_bins.csv', header=header, index=False)
