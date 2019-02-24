# This program plots the detections that make it through into the 
# CDF plots.
import os
import numpy as np
import pandas as pd
from datetime import datetime, date
import matplotlib.pyplot as plt

CATALOG_DIR = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
                'coincident_microbursts_catalogues/')
AC6_DATA_PATH = lambda sc_id: ('/home/mike/research/ac6/ac6{}/'
                                'ascii/level2'.format(sc_id))

class PlotMicrobursts:
    def __init__(self, catalog_version, plot_width=5):
        """
        This class plots the detections from the coincident
        microburst catalog with a set of default filters.
        The user can supply other filter values in the 
        filterDict.
        """
        self.load_catalog(catalog_version)
        return

    def filter_catalog(self, defaultFilter=True, filterDict={}):
        if defaultFilter:
            # High CC filter
            self.catalog = self.catalog[self.catalog['time_cc'] >= 0.8]
            # Rad belt filter
            self.catalog = self.catalog[
                                        (np.abs(self.catalog['Lm_OPQ']) > 4) &
                                        (np.abs(self.catalog['Lm_OPQ']) < 8)
                                        ]
            # Curtain filter
            self.catalog = self.catalog[self.catalog['time_cc'] > 
                                        self.catalog['space_cc']+0.3]
            # USA filter
            self.catalog = self.catalog[
                    ((self.catalog['lon'] > -60) | (self.catalog['lon'] < -140))|
                    ((self.catalog['lat'] > 70) | (self.catalog['lat'] < 15))
                                        ]
            # SAA filter
            self.catalog = self.catalog[
                    ((self.catalog['lon'] > 30)  | (self.catalog['lon'] < -116)) |
                    ((self.catalog['lat'] < -90) | (self.catalog['lat'] > 0))
                                ]
            for key, vals in filterDict.items():
                self.catalog = self.catalog[
                    ((self.catalog[key] > min(vals)) & 
                    (self.catalog[key] < max(vals)))
                                ]
        return

    def load_catalog(self, catalog_version):
        """
        Load the csv microburst catalog. Any key that has the word time in
        it is converted to am array of datetime objects. Time keys are:
        dateTime, time_spatial_A, and time_spatial_B.
        """
        CATALOG_PATH = os.path.join(CATALOG_DIR,
                'AC6_coincident_microbursts_v{}.txt'.format(catalog_version))
        self.catalog = pd.read_csv(CATALOG_PATH)
        # Convert the catalog times to datetime objects
        for timeKey in ['dateTime', 'time_spatial_A', 'time_spatial_B']:
            self.catalog[timeKey] = pd.to_datetime(self.catalog[timeKey])
        return

    def loop(self):
        """
        Loops over the detections that made it through the filters and 
        make space-time plots.
        """
        current_date = date.min

        for _, row in self.catalog.iterrows():
            if row.dateTime.date() != current_date:
                # Load current day AC-6 data if not loaded already
                self.load_ten_hz_data(row.dateTime.date())
                current_date = row.dateTime.date()
        
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
        print(pathA, pathB)
        self.ac6a_data = pd.read_csv(pathA, na_values='-1e+31')
        self.ac6a_data['dateTime'] = pd.to_datetime(self.ac6a_data[time_keys])
        self.ac6b_data = pd.read_csv(pathB, na_values='-1e+31')
        self.ac6b_data['dateTime'] = pd.to_datetime(self.ac6b_data[time_keys])
        return

if __name__ == '__main__':
    p = PlotMicrobursts(4)
    p.filter_catalog(filterDict={'Dist_Total':[100, 200]})
    p.loop()