# Here I calculate the burst parameter for each of the sorted
# microburst detections

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
from datetime import date, datetime
import dateutil.parser
from matplotlib.dates import date2num

sys.path.insert(0, '/home/mike/research/microburst-detection/burst_parameter')
import burst_parameter

AC6_DATA_PATH = lambda sc_id: ('/home/mike/research/ac6/ac6{}/'
                                'ascii/level2'.format(sc_id))

class Microburst_Burst_Parameter:
    def __init__(self, catalog_dir, catalog_name):
        self.catalog_dir = catalog_dir
        self.catalog_name = catalog_name
        return

    def load_catalog(self):
        """
        Load the catalog to loop over.
        """
        catalog_path = os.path.join(self.catalog_dir, self.catalog_name)
        self.catalog = pd.read_csv(catalog_path)
        # Convert the catalog times to datetime objects
        for timeKey in ['dateTime', 'time_spatial_A', 'time_spatial_B']:
            self.catalog[timeKey] = pd.to_datetime(self.catalog[timeKey])

    def loop(self, bp_thresh=3, **burst_param_kwargs):
        """
        Loops over the detections that made it through the filters and 
        make space-time plots.
        """
        current_date = date.min

        self.bp = np.nan*np.ones((self.catalog.shape[0], 2))

        for i, row in self.catalog.iterrows():
            if row.dateTime.date() != current_date:
                # Load current day AC-6 data if not loaded already
                self.load_ten_hz_data(row.dateTime.date())
                current_date = row.dateTime.date()
                # Calculate the burst paramater for that day
                day_bp_A = burst_parameter.obrien_burst_param(
                            self.ac6a_data.dos1rate, 0.1, **burst_param_kwargs)
                day_bp_B = burst_parameter.obrien_burst_param(
                            self.ac6b_data.dos1rate, 0.1, **burst_param_kwargs)

            t0 = date2num(row.dateTime)
            time_an = date2num(self.ac6a_data.dateTime)
            time_bn = date2num(self.ac6b_data.dateTime)
            idta = np.where(time_an >= t0)[0]
            idtb = np.where(time_bn >= t0)[0]
            if bp_thresh:
                self.bp[i, 0] = max(day_bp_A[idta[0]-bp_thresh:idta[0]+bp_thresh])
                self.bp[i, 1] = max(day_bp_B[idtb[0]-bp_thresh:idtb[0]+bp_thresh])
            else:
                self.bp[i, 0] = day_bp_A[idta[0]]
                self.bp[i, 1] = day_bp_B[idtb[0]]
        return

    def save_catalog(self):
        """ Overwrites the input catalog with the burst paramater appended. """
        self.catalog['bp_max'] = self.bp.max(axis=1)
        catalog_path = os.path.join(self.catalog_dir, self.catalog_name)
        self.catalog.to_csv(catalog_path, index=False)
        return

    def load_ten_hz_data(self, day):
        """
        Load the 10 Hz AC6 data from both spacecraft on date.
        """
        time_keys = ['year', 'month', 'day', 'hour', 'minute', 'second']
        dayStr = '{0:%Y%m%d}'.format(day)
        pathA = os.path.join(AC6_DATA_PATH('a'), 
                'AC6-A_{}_L2_10Hz_V03.csv'.format(dayStr))
        pathB = os.path.join(AC6_DATA_PATH('b'), 
                'AC6-B_{}_L2_10Hz_V03.csv'.format(dayStr))
        self.ac6a_data = pd.read_csv(pathA, na_values='-1e+31')
        self.ac6a_data['dateTime'] = pd.to_datetime(self.ac6a_data[time_keys])
        self.ac6b_data = pd.read_csv(pathB, na_values='-1e+31')
        self.ac6b_data['dateTime'] = pd.to_datetime(self.ac6b_data[time_keys])
        return

if __name__ == '__main__':
    catalog_dir = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
                'coincident_microbursts_catalogues')
    catalog_name = 'AC6_coincident_microbursts_sorted_Brady_v6.txt'
    catalog_path = os.path.join(catalog_dir, catalog_name)

    m = Microburst_Burst_Parameter(catalog_dir, catalog_name)
    m.load_catalog()
    m.loop()
    m.save_catalog()

    df = pd.read_csv(catalog_path)
    plt.hist(df.max(axis=1), bins=np.arange(0, 40))
    plt.show()