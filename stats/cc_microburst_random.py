# Cross-correlate microbursts vs random times in the same bin.

# This script finds the 10 most common bins for the microburst CDF as a function of separation.
import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser
import itertools
import os

import pandas as pd

# Microburst catalog path.
version = 5
catalog_path = ('/home/mike/research/ac6_microburst_scale_sizes/'
                'data/microburst_catalogues'
                '/AC6A_microbursts_v{}.txt'.format(version))
# Binned counts directory
bin_counts_dir = ('/home/mike/research/ac6_microburst_scale_sizes/'
            'data/binned_counts')

# Load and filter the microbursts from AC-6 A.
catalog = pd.read_csv(catalog_path)
catalog = catalog[(catalog['Lm_OPQ'] > 4) & (catalog['Lm_OPQ'] < 8)]
# USA filter
catalog = catalog[
    ((catalog['lon'] > -60) | (catalog['lon'] < -140)) |
    ((catalog['lat'] > 70) | (catalog['lat'] < 15))
    ]
catalog['dateTime'] = pd.to_datetime(catalog['dateTime'])

# Load the most common bins file
common_bins = pd.read_csv('most_common_l_mlt_ae_bins.csv')

# Bin widths and helper function to get bin string names.
dL = 1
dMLT = 1
dAE = 100
#N = 20
F = np.nan*np.ones(common_bins.shape[0])

def index_to_name(L, MLT, AE):
    s = '{}_L_{}_{}_MLT_{}_{}_AE_{}'.format(
        int(L), int(L)+dL, int(MLT), int(MLT)+dMLT, int(AE), int(AE)+dAE)
    return s

def filter_catalog(catalog, row):
    """ 
    This function filters the catalog by the common
    bin row L, MLT, and AE values.
    """
    filtered_catalog = catalog[
                            (catalog.Lm_OPQ > row.L) & 
                            (catalog.Lm_OPQ < row.L + dL)
                            ]
    filtered_catalog = filtered_catalog[
                            (filtered_catalog.MLT_OPQ > row.MLT) & 
                            (filtered_catalog.MLT_OPQ < row.MLT + dMLT)
                            ]
    filtered_catalog = filtered_catalog[
                            (filtered_catalog.AE > row.AE) & 
                            (filtered_catalog.AE < row.AE + dAE)
                            ]
    return filtered_catalog

def CC_microburst_random(microbursts, counts, N_CC=1000, N_MAX=10000, 
                        window=10, window_thresh=1, CC_thresh=0.8):
    """ 
    This function takes the microbursts data frame and chooses N_CC 
    number of microbursts to cross correlate (CC) with random counts
    in counts. If data is sparse and N_MAX is reached, the loop will
    abort. The kwargs window and window_thresh define the CC windows. 
    If there is no continious CC window, another random count rate
    time will be picked.
    """
    n_ccd = 0
    n_attempts = 0

    #random_bursts = microbursts.sample(n=N_CC, replace=True)
    CC_arr = np.nan*np.ones(N_CC)

    while n_ccd < N_CC and n_attempts < N_MAX:
        n_attempts += 1
        iA = np.where(np.random.choice(microbursts['dateTime']) == counts['dateTime'])[0]
        if len(iA) == 0:
            continue
        # assert len(iA) == 1, 'None or multiple microburst times found!\n{}\n{}'.format(
        #                         iA, random_bursts.iloc[n_ccd].dateTime)

        iB = np.random.choice(counts.index)

        cc = CC(iA[0], iB, counts, window, window_thresh)
        if not np.isnan(cc):
            CC_arr[n_ccd] = cc
            n_ccd += 1
    signif_cc = len(np.where((CC_arr > CC_thresh) & ~np.isnan(CC_arr))[0])
    if signif_cc == 0:
        return np.nan
    #     print('\nCC_arr =',CC_arr)
    return signif_cc/n_ccd

def CC(iA, iB, counts, CC_width, CC_time_thresh):
    """ 
    This method calculates the normalized cross-correlation 
    between two AC6 time series indexed by iA and iB.
    """
    # Calculate the indicies to cross correlate over
    aStart = iA - CC_width//2 - CC_time_thresh
    aEnd = iA + CC_width//2 + CC_time_thresh
    bStart = iB - CC_width//2
    bEnd = iB + CC_width//2

    # If start or end indicies are outisde of the 
    # self.countsArr index range.
    if aStart < 0 or bStart < 0: 
        return np.nan
    if (aEnd >= len(counts['dateTime']) or 
            bEnd >= len(counts['dateTime'])): 
        return np.nan

    # If start and end indicies are not "close" to each
    # other, return error code np.nan.
    dtA = round((counts['dateTime'].loc[aEnd] - 
                counts['dateTime'].loc[aStart]).total_seconds(), 1)
    dtB = round((counts['dateTime'].loc[bEnd] - 
                counts['dateTime'].loc[bStart]).total_seconds(), 1)
    if dtA > (CC_width+2*CC_time_thresh)/10: return np.nan
    if dtB > (CC_width)/10: return np.nan

    # Cross-correlate starting here
    norm = np.sqrt(len(counts.loc[aStart:aEnd, 'dos1rate'])*\
                    len(counts.loc[bStart:bEnd, 'dos1rate'])*\
                    np.var(counts.loc[aStart:aEnd, 'dos1rate'])*\
                    np.var(counts.loc[bStart:bEnd, 'dos1rate'])) 
    # Mean subtraction.
    x = (counts.loc[aStart:aEnd, 'dos1rate'] - 
        counts.loc[aStart:aEnd, 'dos1rate'].mean() )
    y = (counts.loc[bStart:bEnd, 'dos1rate'] - 
        counts.loc[bStart:bEnd, 'dos1rate'].mean() )
    # Cross-correlate
    ccArr = np.correlate(x, y, mode='valid')
    # Normalization
    ccArr /= norm
    return max(ccArr)

# Loop over the N most common bins
for index, row in common_bins.iterrows():
# for index, row in common_bins.iloc[0:10, :].iterrows():
    # Find all of the microbursts in that bin
    filtered_catalog = filter_catalog(catalog, row)
    
    # Load the count bin
    s = index_to_name(row.L, row.MLT, row.AE)
    binName = 'AC6_counts_' + s + '.csv'
    binPath = os.path.join(bin_counts_dir, binName)
    binDf = pd.read_csv(binPath, names=['dateTime', 'dos1rate'])
    binDf['dateTime'] = pd.to_datetime(binDf['dateTime'])

    F[index] = CC_microburst_random(filtered_catalog, binDf, 
                                    N_CC=1000, N_MAX=10000, 
                                    window=10, window_thresh=1, 
                                    CC_thresh=0.8)
print(F)
print(100*np.sqrt(np.nansum(F**2)))