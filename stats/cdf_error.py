import collections
import os

import pandas as pd
class CDF_error:
    def __init__(self, catV, common_bin_path=None, 
                 catalog_path=None, bin_dir=None, bin_name=None):
        """
        This class calculates the CDF error from L-MLT-AE bins
        defined in the most_common_count_bins.csv file.
        """
        # Determine appropriate paths
        # Determine the most common bin paths. 
        if common_bin_path is None:
            bin_dir = '/home/mike/research/ac6_microburst_scale_sizes/data'
            bin_name = 'most_common_counts_bins.csv'
            common_bin_path = os.path.join(bin_dir, bin_name)
        # Determine the coincident microburst catalog path.
        if catalog_path is None:
            catalog_dir = ('/home/mike/research/ac6_microburst_scale_sizes/'
                            'data/coincident_microbursts_catalogues')
            catalog_name = f'AC6_coincident_microbursts_sorted_v{catV}.txt'
            catalog_path = os.path.join(catalog_dir, catalog_name)
        # Determine the counts bin file directory and naming format.
        self.dL = 1
        self.dMLT = 1
        self.dAE = 100
        if bin_dir is None:
            self.bin_dir = ('/home/mike/research/ac6_microburst_scale_sizes/'
                            'data/binned_counts')
        else:
            self.bin_dir = bin_dir

        if bin_name is None:
            self.bin_name = lambda L, MLT, AE: (f'AC6_counts_{int(L)}_L_{int(L+self.dL)}'
                                                f'_{int(MLT)}_MLT_{int(MLT+self.dMLT)}'
                                                f'_{int(AE)}_AE_{int(AE+self.dAE)}.csv')
        else:
            self.bin_name = bin_name
            
        # Load the most common bins
        self.load_common_count_bins(common_bin_path)
        # Load the catalog
        self.load_microburst_catalog(catalog_path)
        return

    def load_common_count_bins(self, path):
        """ Load and process most common count bin csv file. """
        self.count_bins = collections.defaultdict(list)
        # Read csv file.
        df = pd.read_csv(path)
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

if __name__ == '__main__':
    err = CDF_error(6)