import csv
import numpy as np
import itertools
import dateutil.parser
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import os
import sys

import mission_tools.ac6.read_ac_data as read_ac_data

class PlotCrossCorrelatedBursts:
    def __init__(self, catPath, saveDir='/home/mike/Desktop/validation_plots'):
        self.saveDir = os.path.join(saveDir, datetime.now().date().isoformat())
        if not os.path.exists(self.saveDir):
            os.mkdir(self.saveDir)
            print('Made directory:', self.saveDir)

        self.cat = self._load_catalog_(catPath)
        return

    def filter_catalog(self, flt):
        """ 
        This method filters self.cat according to the filters in flt.
        Filter is a tuple of length 3 with elements (key, lower_bound, 
        upper_bound). The key is the data key to filter by, and lower_bound
        and upper_bound values are for the optional upper and lower bound
        values to filter for that key.
        """
        if flt[1] is None:
            idx = np.where(self.cat[flt[0]] < flt[2])[0]
        elif flt[2] is None:
            idx = np.where(self.cat[flt[0]] > flt[1])[0]
        else:
            idx = np.where(
                            (self.cat[flt[0]] > flt[1]) & 
                            (self.cat[flt[0]] < flt[2])
                            )[0]
        for key in self.cat:
            self.cat[key] = self.cat[key][idx]
        return

    def loop(self, pltWidth=5, ccWidth=1):
        """ 
        This method will loop over the catalog and make plots of 
        coincident microbursts at the same time and position. CC
        information will be annotated.
        """
        prev_date = datetime.min
        _, self.ax = plt.subplots(2, figsize=(8, 8))

        for row, t0 in enumerate(self.cat['dateTime']):
            # Load data on t0 if it has not been loaded already.
            if t0.date() != prev_date.date():
                print('Loading data from', t0.date())
                self.dataA = read_ac_data.read_ac_data_wrapper('A', t0)
                self.dataB = read_ac_data.read_ac_data_wrapper('B', t0)
                prev_date = t0

            # Plot count rates.
            plotdt = timedelta(seconds=pltWidth/2)
            ccdt = timedelta(seconds=ccWidth/2)

            idA = np.where( (self.dataA['dateTime'] > t0-plotdt) & 
                            (self.dataA['dateTime'] < t0+plotdt) &
                            (self.dataA['dos1rate'] >= 0)
                            )[0]
            idB = np.where( (self.dataB['dateTime'] > t0-plotdt) & 
                            (self.dataB['dateTime'] < t0+plotdt) &
                            (self.dataB['dos1rate'] >= 0)
                            )[0] 
            cA = self.dataA['dos1rate'][idA] - self.dataA['dos1rate'][idA].mean()
            cB = self.dataB['dos1rate'][idB] - self.dataB['dos1rate'][idB].mean()

            self.ax[0].plot(self.dataA['dateTime'][idA], cA, 'r', 
                    label='AC6A dos1rate')
            self.ax[0].plot(self.dataB['dateTime'][idB], cB, 'b', 
                    label='AC6B dos1rate')
            self.ax[0].axvline(t0, c='k')
            self.ax[0].text(0.05, 0.9, 'time_CC = {}'.format(round(self.cat['time_cc'][row], 2)), 
                            transform=self.ax[0].transAxes, ha='left')
            self.ax[0].legend(loc=1)
            self.ax[0].set_ylabel('Counts/s')
            self.ax[1].set_xlabel('UTC')

            # Save plot.
            saveDate = t0.replace(microsecond=0).isoformat().replace(':', '')
            saveName = '{}_microburst_validation.png'.format(saveDate)
            plt.savefig(os.path.join(self.saveDir, saveName))

            for a in self.ax:
                a.clear() # Clear axis.

        return

    def _load_catalog_(self, fPath):
        """
        This method reads in a catalog csv file, and saves it to a dictionary.
        """
        dtypes = [datetime] + 14*[float] + 2*[datetime]
        converters = {0:dateutil.parser.parse,
                -2:dateutil.parser.parse,
                -1:dateutil.parser.parse}

        rawData = np.genfromtxt(fPath, delimiter=',', names=True, dtype=dtypes,  
                            converters=converters)
        data = {}
        for col, key in enumerate(rawData.dtype.names):
            data[key] = np.array([row[col] for row in rawData])
        return data 

if __name__ == '__main__':
    path = './../stats/coincident_microburst_test_v2.csv'
    p = PlotCrossCorrelatedBursts(path)
    p.filter_catalog(('Dist_Total', 60, 100))
    p.filter_catalog(('time_cc', 0.79, None))
    p.loop()