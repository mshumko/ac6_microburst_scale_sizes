# This program plots the detections that make it through into the 
# CDF plots.
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CATALOG_DIR = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
                'coincident_microbursts_catalogues/')

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
        for timeKey in filter(lambda key: 'time' in key.lower(), self.catalog.keys()):
            self.catalog[timeKey] = pd.to_datetime(self.catalog[timeKey])
        return

    def loop(self):

        return

    def load_ten_hz_data(self, date):
        """
        Load the 10 Hz AC6 data from both spacecraft on date.
        """

        return

if __name__ == '__main__':
    p = PlotMicrobursts(4)