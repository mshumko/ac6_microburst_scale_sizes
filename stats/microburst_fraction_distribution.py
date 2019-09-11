import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import dateutil.parser
import os

import mission_tools.ac6.read_ac_data as read_ac_data
import IRBEM

Re = 6371 # Earth radius, km

class LeoMicroburstFraction:
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
            print(f'Progress {round(i/(self.bins.shape[0]-1), 2)}%') 
            # Filter the full catalog in separation and L shell.
            filtered_catalog = self.microburst_catalog[
                (self.microburst_catalog.Dist_Total > bin_i) &
                (self.microburst_catalog.Dist_Total < bin_f) &
                (self.microburst_catalog.Lm_OPQ > 4) &
                (self.microburst_catalog.Lm_OPQ < 8) &
                # There is a valid temporal CC value (data exists for both spacecraft)
                (~np.isnan(self.microburst_catalog.time_cc)) 
                ]
            # Now find the number of detections made by one AC6 unit and NOT
            # the other.
            n[i] = filtered_catalog.shape[0]
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

    def plot_fraction(self):
        """ Plot fraction of microbursts with steps and vertical bars for error bars."""
        plt.step(self.bins[:-1], self.f, where='post')
        plt.errorbar(self.bins[:-2]+self.bin_width/2, self.f[:-1], ls='', yerr=self.f_err[:-1])
        plt.xlabel('AC6 separation [km]')
        plt.ylabel('Fraction of microbursts')
        plt.title(f'AC6 Coincident Microburst Probability\n 4 < L < 8')
        return

class EquatorialMicroburstFraction(LeoMicroburstFraction):
    def __init__(self, sc_id, microburst_catalog_name, c_microburst_catalog_name,
                microburst_catalog_dir=None, c_microburst_catalog_dir=None):
        """ 
        Child class of LeoMicroburstFraction except the make_microburst_fraction 
        and plot_fraction methods are overwritten and map to the magnetic 
        equator using OPQ magnetic field model.
        """
        super().__init__(sc_id, microburst_catalog_name, c_microburst_catalog_name,
                        microburst_catalog_dir=microburst_catalog_dir, 
                        c_microburst_catalog_dir=c_microburst_catalog_dir)
        self.model = IRBEM.MagFields(kext='OPQ77')
        return

    def make_microburst_fraction(self, key='d_equator',
                                start_bin=0, end_bin=1000, bin_width=100):
        """ 
        Overwritten class where all of the detections are first mapped 
        to the magnetic equator.
        """
        # First map all events to the magnetic equator.
        if not 'd_equator' in self.microburst_catalog.columns:
            self.microburst_catalog.loc[:, 'd_equator'] = np.array([
                        self._map2equator(row.lat, row.lon, row.alt, 
                                        row.dateTime, row.Dist_Total) 
                        for _, row in self.microburst_catalog.iterrows()])
        # Then map the coincident events to the magnetic equator.
        if not 'd_equator' in self.c_microburst_catalog.columns:
            self.c_microburst_catalog.loc[:, 'd_equator'] = np.array([
                        self._map2equator(row.lat, row.lon, row.alt, 
                                        row.dateTime, row.Dist_Total) 
                        for _, row in self.c_microburst_catalog.iterrows()])
        # Calculate fraction
        self.bin_width = bin_width
        self.bins = np.arange(start_bin, end_bin+1, bin_width)
        self._microburst_fraction(key=key)
        return

    def plot_fraction(self):
        """ Plot fraction of microbursts with steps and vertical bars for error bars."""
        plt.step(self.bins[:-1], self.f, where='post')
        plt.errorbar(self.bins[:-2]+self.bin_width/2, self.f[:-1], ls='', yerr=self.f_err[:-1])
        plt.xlabel('AC6 equatorial separation [km]')
        plt.ylabel('Fraction of microbursts')
        plt.title(f'Equatorial AC6 Coincident Microburst Probability\n 4 < L < 8')
        plt.xlim(self.bins[0], self.bins[-1])
        return

    def save_fraction(self, file_path=None):
        """ Saves the microburst fraction and errors to csv file. """
        if file_path is None:
            file_dir = '/home/mike/research/ac6_microburst_scale_sizes/data'
            file_name = 'equatorial_microburst_fraction.csv'
            file_path = os.path.join(file_dir, file_name)
        df = pd.DataFrame({'f':self.f, 'f_err':self.f_err})
        df.to_csv(file_path, index=False)
        return

    def _map2equator(self, lat, lon, alt, time, d):
        """ Maps to magnetic equator assuming d is soly in latitude. """
        # Define the coordinates of the two spacecraft.
        dLat = self._deltaLat(d, alt)
        X1 = {'x1':alt, 'x2':lat-dLat, 'x3':lon, 'dateTime':time}
        X2 = {'x1':alt, 'x2':lat+dLat, 'x3':lon, 'dateTime':time}
        # Run IRBEM
        X1_equator = self.model.find_magequator(X1, None)['XGEO']
        X2_equator = self.model.find_magequator(X2, None)['XGEO']
        # Calculate the separations.
        self.d_equator = Re*np.linalg.norm(X1_equator-X2_equator)
        return self.d_equator

    def _deltaLat(self, d, alt):
        """
        Calculate the half of the change in angle for a spacecraft at
        an altitude alt and separated by a distance d.
        """
        return np.rad2deg(np.arcsin(d/(2*(Re+alt))))

    def _microburst_fraction(self, key='d_equator'):
        """ 
        This method calculates the fraction of coincident to all 
        microbursts in separation bins defined by the bin attribute. The
        separations are estimated using the key kwarg
        """
        n = np.nan*np.ones(len(self.bins)-1)
        n_c = np.nan*np.ones(len(self.bins)-1)

        # Loop over the separation bins and calculate the number of total 
        # and coincident microburst events.
        for i, (bin_i, bin_f) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # Filter the full catalog in separation and L shell.
            filtered_catalog = self.microburst_catalog[
                (self.microburst_catalog[key] > bin_i) &
                (self.microburst_catalog[key] < bin_f) &
                (self.microburst_catalog.Lm_OPQ > 4) &
                (self.microburst_catalog.Lm_OPQ < 8) &
                # There is a valid temporal CC value (data exists for both spacecraft)
                (~np.isnan(self.microburst_catalog.time_cc)) 
                ]

            # Now find the number of microbursts observed by AC6B only.
            n = 0
            for _, row in filtered_catalog.iterrows():
                if row.dateTime == row.time_spatial_B:
                    n += 1

            #n[i] = filtered_catalog.shape[0]
            n_c[i] = np.sum(
                        (self.c_microburst_catalog[key] > bin_i) &
                        (self.c_microburst_catalog[key] < bin_f)
                        )
        self.f = n_c/n
        # Need to check how the uncertanity is calculated.
        # Fraction error propagation eq. from wiki 
        # https://en.wikipedia.org/wiki/Propagation_of_uncertainty
        self.f_err = self.f*np.sqrt( 1/n_c + 1/n )
        return self.f, self.f_err
 


if __name__ == '__main__':
    sc_id = 'AC6B'
    microburst_name = 'AC6_coincident_microbursts_v8.txt' #'AC6B_microbursts_v5.txt'
    microburst_catalog_dir = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
                              'coincident_microbursts_catalogues')
    # coincident_catalog_name = 'AC6_coincident_microbursts_sorted_Brady_v6.txt'
    coincident_catalog_name = 'AC6_coincident_microbursts_sorted_v6.txt'

    mf = EquatorialMicroburstFraction(sc_id, microburst_name, coincident_catalog_name, 
                                    microburst_catalog_dir=microburst_catalog_dir)
    mf.make_microburst_fraction(start_bin=0, end_bin=2000, bin_width=100)
    mf.save_fraction()
    mf.plot_fraction()
    plt.show()