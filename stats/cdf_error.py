import collections
import os
import numpy as np

import pandas as pd

class CDF_error:
    def __init__(self, catV=None, common_bin_path=None, 
                 catalog_path=None, count_bin_dir=None, count_bin_name=None):
        """
        This class calculates the CDF error from L-MLT-AE bins
        defined in the most_common_count_bins.csv file.
        """
        # Determine appropriate paths
        # Determine the most common bin paths. 
        if common_bin_path is None:
            bin_dir = '/home/mike/research/ac6_microburst_scale_sizes/data/'
            bin_name = 'most_common_counts_bins.csv'
            common_bin_path = os.path.join(bin_dir, bin_name)
        # Determine the coincident microburst catalog path.
        if catalog_path is None:
            # catalog_dir = ('/home/mike/research/ac6_microburst_scale_sizes/'
            #                 'data/coincident_microbursts_catalogues')
            # catalog_name = f'AC6_coincident_microbursts_sorted_v{catV}.txt'
            catalog_dir = ('/home/mike/research/ac6_microburst_scale_sizes/'
                            'data/microburst_catalogues')
            catalog_name = f'AC6A_microbursts_v{catV}.txt'
            catalog_path = os.path.join(catalog_dir, catalog_name)
        # Determine the counts bin file directory and naming format.
        self.dL = 1
        self.dMLT = 1
        self.dAE = 100
        if count_bin_dir is None:
            self.count_bin_dir = ('/home/mike/research/ac6_microburst_scale_sizes/'
                            'data/binned_counts')
        else:
            self.count_bin_dir = bin_dir

        if count_bin_name is None:
            self.count_bin_name = lambda L, MLT, AE: (f'AC6_counts_{int(L)}_L_{int(L+self.dL)}' +
                                                    f'_{int(MLT)}_MLT_{int(MLT+self.dMLT)}' +
                                                    f'_{int(AE)}_AE_{int(AE+self.dAE)}.csv')
        else:
            self.count_bin_name = count_bin_name

        # Load the most common bins
        self.load_common_count_bins(common_bin_path)
        # Load the catalog
        self.load_microburst_catalog(catalog_path)
        return

    def loop(self, N_CC=1000, N_MAX=10000, window=10, 
            window_thresh=1, CC_thresh=0.8):
        """ 
        Main loop that does the cross correlating between
        microbursts and random count times
        """
        count_bin_keys = list(self.count_bins)
        # Define an array F that is n_sep_bins x n_common_bins in size.
        self.F = np.nan*np.ones((len(count_bin_keys), 
                                len(self.count_bins[count_bin_keys[0]])))
        # Loop over separation bins
        for i, key_d in enumerate(count_bin_keys):
            # Loop over the most common bins in each separation bin.
            for j, row_i in enumerate(self.count_bins[key_d]):
                # Pass the counts bin and microbursts in that bin to
                # self.CC_wrapper to cross-correlate.
                counts_bin_path = os.path.join(self.count_bin_dir, 
                        self.count_bin_name(row_i.L, row_i.MLT, row_i.AE))
                # counts = pd.read_csv(counts_bin_path)
                counts = pd.read_csv(counts_bin_path, names=['dateTime', 'dos1rate'])
                counts['dateTime'] = pd.to_datetime(counts['dateTime'])

                microbursts = self._get_microbursts(key_d, row_i)

                if microbursts.shape[0] > 100:
                    self.F[i, j] = self.CC_wrapper(microbursts, counts, 
                                    N_CC, N_MAX, window, window_thresh,
                                    CC_thresh)
        return

    def CC_wrapper(self, microbursts, counts, N_CC, N_MAX, window, 
                    window_thresh, CC_thresh):
        """ 
        This method finds the microburst detections made above
        a separation d and L-MLT-AE bin befined by row and cross
        correlates it a bunch of times.
        """
        n_ccd = 0
        n_attempts = 0
        CC_arr = np.nan*np.ones(N_CC)

        while n_ccd < N_CC and n_attempts < N_MAX:
            n_attempts += 1
            iA = np.where(np.random.choice(microbursts['dateTime']) == counts['dateTime'])[0]
            if len(iA) == 0:
                continue

            iB = np.random.choice(counts.index)

            cc = self.CC(iA[0], iB, counts, window, window_thresh)
            if not np.isnan(cc):
                CC_arr[n_ccd] = cc
                n_ccd += 1
                
        signif_cc = len(np.where((CC_arr > CC_thresh) & ~np.isnan(CC_arr))[0])
        if signif_cc == 0:
            return np.nan
        #     print('\nCC_arr =',CC_arr)
        return signif_cc/n_ccd

    def CC(self, iA, iB, counts, CC_width, CC_time_thresh):
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

    def load_common_count_bins(self, path):
        """ Load and process most common count bin csv file. """
        self.count_bins = collections.defaultdict(list)
        # Read csv file.
        df = pd.read_csv(path, dtype={'L':int, 'MLT':int, 'AE':int})
        # Loop over the count bin rows
        for _, row in df.iterrows():
            row_filtered = row.drop('d')
            self.count_bins[row.d].append(row_filtered)
        return

    def load_microburst_catalog(self, path, defaultFilter=True):
        """ 
        Load the microburst catalog file and optionally apply 
        a default filter.
        """
        self.catalog = pd.read_csv(path)
        if defaultFilter:
            # L filter
            self.catalog = self.catalog[(self.catalog['Lm_OPQ'] > 4) & 
                                        (self.catalog['Lm_OPQ'] < 8)]
        # Filter out data over the USA
        self.catalog = self.catalog[
            ((self.catalog['lon'] > -60) | (self.catalog['lon'] < -140)) |
            ((self.catalog['lat'] > 70)  | (self.catalog['lat'] < 15))
            ]
        # Convert times.
        #catalog['dateTime'] = pd.to_datetime(catalog['dateTime'])
        return

    def _get_microbursts(self, d, row):
        """ 
        Helper function to filter microbursts out of self.catalog
        """
        m = self.catalog.copy()
        m = m[m['Dist_Total'] > d]
        m = m[m['Lm_OPQ'] > row.L]
        m = m[m['Lm_OPQ'] < row.L + self.dL]
        m = m[m['MLT_OPQ'] > row.MLT]
        m = m[m['MLT_OPQ'] < row.MLT + self.dMLT]
        m = m[m['AE'] > row.AE]
        m = m[m['AE'] < row.AE + self.dAE]
        return m

if __name__ == '__main__':
    err = CDF_error(catV=5)
    err.loop()