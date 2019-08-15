import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import dateutil.parser
import os

import mission_tools.ac6.read_ac_data as read_ac_data

class MicroburstFraction:
    def __init__(self, sc_id, microburst_catalog_name, c_microburst_catalog_name,
                microburst_catalog_dir=None, c_microburst_catalog_dir=None):
        """
        This class calculates and plots the fraction of simulatenous to 
        all microburst detections as a function of separation. Uncertanity 
        is first assumed due to Poisson noise and will be then expanded to
        a systematic uncertanity due to transmitter noise.
        """
        self.sc_id = sc_id

        # Specify the coincident catalog directory
        if c_microburst_catalog_dir is None:
            self.c_microburst_catalog_dir = ('/home/mike/research/'
                            'ac6_microburst_scale_sizes/'
                            'data/coincident_microbursts_catalogues')
        else:
            self.c_microburst_catalog_dir = c_microburst_catalog_dir

        # Specify the entire microburst catalog directory
        if microburst_catalog_dir is None:
            self.microburst_catalog_dir = ('/home/mike/research/'
                            'ac6_microburst_scale_sizes/data/'
                            'microburst_catalogues')
        else:
            self.microburst_catalog_dir = microburst_catalog_dir
        
        # Load the catalogs
        self.microburst_catalog = self.load_catalog(
            os.path.join(self.microburst_catalog_dir, microburst_catalog_name))
        self.c_microburst_catalog = self.load_catalog(
            os.path.join(self.c_microburst_catalog_dir, c_microburst_catalog_name))
        return

    def load_catalog(self, path):
        """ 
        Loads a microburst catalog (or) any CSV file with column
        names in header. Path is the absolute path to the catalog
        file. Returns a pandas DataFrame object with the catalog. 
        """
        return pd.read_csv(path)

    def make_microburst_fraction(self, systematic_error=None,
                                start_bin=0, end_bin=100, bin_width=5):
        """ 
        This method calculates the fraction of coincident to all 
        microbursts in separation bins defined by start/end_bin and 
        bin_width kwargs. The systematic_error is the error that is added
        to the statistical error in quadrature.
        """
        self.bins = np.arange(start_bin, end_bin+1, bin_width)
        self.bin_width = bin_width
        n = np.nan*np.ones(len(self.bins)-1)
        n_c = np.nan*np.ones(len(self.bins)-1)

        # Loop over the separation bins and calculate the number of total 
        # and coincident microburst events.
        for i, (bin_i, bin_f) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # NEED TO ADD LOOP FOR n[i] to only count microbursts observed by
            # one spacecraft AND the other spacecraft did not see anything 
            # (but had data!) 
            filtered_catalog = self.microburst_catalog[
                (self.microburst_catalog.Dist_Total > bin_i) &
                (self.microburst_catalog.Dist_Total < bin_f) &
                (self.microburst_catalog.Lm_OPQ > 4) &
                (self.microburst_catalog.Lm_OPQ < 8) 
                ]
            # Now find the number of detections made by one AC6 unit and NOT
            # the other.
            n[i] = self.mutually_exclusive_detections(filtered_catalog)
            n_c[i] = np.sum(
                        (self.c_microburst_catalog.Dist_Total > bin_i) &
                        (self.c_microburst_catalog.Dist_Total < bin_f)
                        )
        self.f = n_c/n
        # Need to check how the uncertanity is calculated.
        if systematic_error is None:
            # Fraction error propagation eq. from wiki 
            # https://en.wikipedia.org/wiki/Propagation_of_uncertainty
            self.f_err = self.f*np.sqrt( 1/n_c + 1/n )
        else:
            self.f_err = np.sqrt(
                (self.f*(1-self.f)/np.sqrt(n))**2 + systematic_error**2
                                )
        return

    def mutually_exclusive_detections(self, catalog):
        """ 
        For each detection, check if the other spacecraft had data. Returns
        the number of microburst detections in the catalog AND other sc had
        data at that time.
        """
        n = 0
        current_date = datetime.min
        dt = timedelta(seconds=0.5)
        # Determine the other spacecraft
        if self.sc_id.upper() == 'AC6A':
            other_sc_id = 'B'
        else:
            other_sc_id = 'A'

        for _, row in catalog.iterrows():
            t0 = dateutil.parser.parse(row.dateTime)
            # If the current date is not loaded, load it in and set current_date
            if t0.date != current_date:
                try:
                    self.ac6_10hz = read_ac_data.read_ac_data_wrapper(other_sc_id, 
                                    t0, dType='10Hz')
                    current_date = t0.date()
                
                # Move on if there is something wrong with the file (does not 
                # exist, or empty)
                except AssertionError as err:
                    if str(err) == 'File is empty!':
                        continue
                    else:
                        raise

            ac6_data_flt = self.ac6_10hz[
                                        (self.ac6_10hz.dateTime > t0-dt) & 
                                        (self.ac6_10hz.dateTime < t0+dt)
                                        ]
            if ac6_data_flt.shape[0]:
                n += 1
        return n

    def plot_fraction(self):
        """ Plot fraction of microbursts with steps and vertical bars for error bars."""
        plt.step(mf.bins[:-1], mf.f, where='post')
        plt.errorbar(mf.bins[:-1]+self.bin_width/2, mf.f, ls='', yerr=self.f_err)
        return

if __name__ == '__main__':
    sc_id = 'AC6B'
    microburst_name = 'AC6B_microbursts_v5.txt'
    # coincident_catalog_name = 'AC6_coincident_microbursts_sorted_Brady_v6.txt'
    coincident_catalog_name = 'AC6_coincident_microbursts_sorted_v6.txt'

    mf = MicroburstFraction(sc_id, microburst_name, coincident_catalog_name)
    mf.make_microburst_fraction()
    mf.plot_fraction()
    plt.show()