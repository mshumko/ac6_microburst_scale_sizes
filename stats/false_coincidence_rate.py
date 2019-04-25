# This program calculates the false microburst rate for coincident 
# microbursts.

import numpy as np
import os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

import pandas as pd

class False_Coincidence:
    def __init__(self, coincident_name, microburst_a_name, microburst_b_name, save_path):
        """
        Given a dataset from input_path, calculate the random 
        microburst coicidence rate for each coincident microburst
        identified in input_path.
        """
        coincident_dir = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
                        'coincident_microbursts_catalogues')
        microburst_dir = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
                        'microburst_catalogues')
        self.coincident_microbursts = self._load_catalog(
                        os.path.join(coincident_dir, coincident_name))
        self.microbursts_a = self._load_catalog(
                        os.path.join(microburst_dir, microburst_a_name))
        self.microbursts_b= self._load_catalog(
                        os.path.join(microburst_dir, microburst_b_name))
        self.save_path = save_path
        return

    def loop(self, microburst_width_seconds=0.5, time_intergration_seconds=60):
        """ 
        Loops through the coincident microbursts and calculates the 
        mean occurance rate of microbursts in time_intergration_seconds time 
        span around that microburst time.
        """
        self.coincident_microbursts['false_rate'] = np.nan
        # Probably a better way to do this...
        for index in self.coincident_microbursts.index:
            self.coincident_microbursts['false_rate'].loc[index] = self._get_rate(
                index, microburst_width_seconds, time_intergration_seconds)
        return

    def save_catalog(self):
        """ Saves the coincident microburst dataframe to a csv file. """
        self.coincident_microbursts.to_csv(self.save_path)
        return

    def _load_catalog(self, path):
        """ Loads the microburst catalog """
        catalog = pd.read_csv(path, index_col=0)
        catalog.index = pd.to_datetime(catalog.index)
        return catalog

    def _get_rate(self, t0, microburst_width_seconds, time_intergration_seconds):
        """ 
        Helper function to calculate the mean microburst occurance centered 
        on time t0 and full range time_range_seconds. Microbursts are assumed
        to all be the same width of microburst_width_seconds
        """
        dt = pd.Timedelta(seconds=time_intergration_seconds/2)
        # Find all microbursts in AC6-A and AC6-B data in the time range.
        bursts_a = self.microbursts_a.loc[t0-dt:t0+dt]
        bursts_b = self.microbursts_b.loc[t0-dt:t0+dt]
        # Take the mean number of microbursts and calculate the microburst rate.
        # mean_number_bursts_interval = (bursts_a.shape[0] + bursts_b.shape[0])/2
        # rate = microburst_width_seconds*mean_number_bursts_interval/time_intergration_seconds
        weight = (microburst_width_seconds/time_intergration_seconds)**2
        false_rate = bursts_a.shape[0]*bursts_b.shape[0]*weight
        return false_rate

if __name__ == '__main__':
    coincident_name = 'AC6_coincident_microbursts_sorted_v6.txt'
    microburst_a_name = 'AC6A_microbursts_v5.txt'
    microburst_b_name = 'AC6B_microbursts_v5.txt'

    save_dir = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
            'coincident_microbursts_catalogues')
    save_name = 'AC6_coincident_microbursts_sorted_err_v6.txt'
    save_path = os.path.join(save_dir, save_name)

    f = False_Coincidence(coincident_name, microburst_a_name, 
                            microburst_b_name, save_path)
    f.loop()
    f.save_catalog()

    plt.hist(f.coincident_microbursts.false_rate*100, bins=np.linspace(0, 15, num=25))
    plt.title('AC6 microburst false coincidence rate\n'
             'microburst_width = 0.5 s | integration_width = 60 s')
    plt.xlabel('Rate [%]')
    plt.ylabel('Number of Microbursts')
    plt.show()