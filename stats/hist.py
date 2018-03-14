# This function calculates the normalization parameters for the 
# micorburst scale size distributions.

from datetime import datetime, timedelta
import csv
import numpy as np
import sys

# Import personal libraries
sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

class Hist1D:
    def __init__(self, d=None, sc_id='B', startDate=datetime(2014, 1, 1),
                 endDate=datetime.now(), filterDict={}):
        """
        This class calculates the 1D histograms as a function of distance
        for the various filter parameters. 

        d is the distance bin edges
        """
        if d is None:
            self.d = np.arange(0, 501, 10)
        else:
            self.d = d
        self.count = np.zeros(len(self.d)-1) # Number of days at that separation.
        dDays = (endDate - startDate).days
        self.dates = [startDate + timedelta(days=i) for i in range(dDays)] 
        self.sc_id = sc_id  
        self.filterDict = filterDict
        return

    def __iter__(self):
        """
        This class will iterate over the dats from startDate to endDate, 
        load data, filter it, and histgram it.
        """
        return self

    def get_data(self, date):
        """
        This generator function will load in one day of data at a 
        time, and return control to user.
        """
        try:   
            print('Loading data from {}'.format(date))
            self.ac6data = read_ac_data.read_ac_data_wrapper(
                self.sc_id, date, dType='10Hz', plot=False)
        except AssertionError as err:
            if ( ('None or > 1 AC6 files found' in str(err)) or
                ('Error, the data is not 2D (Empty file?)' in str(err)) ):
                return None
            else:
                raise
        return self.ac6data

    def filterData(self):
        """
        This function will filter the AC-6 data.
        """
        ind = len(self.ac6data[self.ac6data.keys()[0]])
        for key, value in self.filterDict.items():
            idx = np.where((self.ac6data[key] > value.min()) & 
                            (self.ac6data[key] < value.max()))[0]
            ind = np.intersect1d(ind, idx)
        return ind


class Hist2D(Hist1D):
    def __init__(self):
        """
        This class calculates the 2D histograms as a function of distance
        for the various filter parameters. 
        """
        pass

if __name__ == '__main__':
    h = Hist1D()
    #print(f)
    for i, d in enumerate(h):
        print(i)
        print(d)