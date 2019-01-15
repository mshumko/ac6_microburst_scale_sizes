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
        self.cat = self._load_catalog_(catPath)
        
        saveDir = sys.path.join( (saveDir, datetime.now().date()) )
        if not os.path.exists(saveDir):
            os.mkdir(saveDir)
            print('Made directory:', saveDir)
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

    def loop(self):
        """ 
        This method will loop over the catalog and make plots of 
        coincident microbursts at the same time and position. CC
        information will be annotated.
        """

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