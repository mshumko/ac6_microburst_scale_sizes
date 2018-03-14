# This function calculates the normalization parameters for the 
# micorburst scale size distributions.

from datetime import datetime, timedelta
import csv
import numpy as np
import sys
import os

# Import personal libraries
sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

class Hist1D:
    def __init__(self, d=None, sc_id='B', startDate=datetime(2014, 1, 1),
                 endDate=datetime.now(), filterDict={}, flag=True):
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
        self.flag = flag
        return

    def loop_data(self):
        """

        """
        for day in self.dates:
            self.load_day_data(day)
            if self.ac6data is None:
                continue # If data file is empty
            ind = self.filterData()
            self.hist_data(ind)

        return

    def load_day_data(self, date):
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
                ('File is empty!' in str(err)) ):
                self.ac6data = None
                return None
            else:
                raise
        return self.ac6data

    def filterData(self):
        """
        This function will filter the AC-6 data.
        """
        if self.flag:
            ind = np.where(self.ac6data['flag'] == 0)[0]
        else:
            ind = range(len(self.ac6data['flag']))

        for key, value in self.filterDict.items():
            idx = np.where((self.ac6data[key] > np.min(value)) & 
                            (self.ac6data[key] < np.max(value)))[0]
            ind = np.intersect1d(ind, idx)
        return ind

    def hist_data(self, ind):
        """
        This method will histrogram the total distance data.
        """
        H, _ = np.histogram(self.ac6data['Dist_Total'][ind], bins=self.d)
        self.count += H/10

        print(self.count)
        return

    def save_data(self, fPath):
        """
        This method saves the normalization data to a csv file.
        """
        with open(fPath, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['Separation [km]', 'Seconds'])

            for (d, s) in zip(self.d, self.count):
                w.writerow([d, s])
        print('Saved data to {}'.format(os.path.basename(fPath)))
        return

class Hist2D(Hist1D):
    def __init__(self):
        """
        This class calculates the 2D histograms as a function of distance
        for the various filter parameters. 
        """
        pass

if __name__ == '__main__':
    # This script will loop over the intervals in Lbins, and 
    # make normalization histograms for all of them.
    Lbins = [3, 4, 5, 6, 7]
    for (lL, uL) in zip(Lbins[:-1], Lbins[1:]):
        h = Hist1D(filterDict={'Lm_OPQ':[lL,uL]})
        h.loop_data()
        h.save_data('/home/mike/research/ac6-microburst-scale-sizes'
                    '/data/norm/ac6_norm_{}_L_{}.txt'.format(lL, uL))

    # h = Hist1D()
    # h.loop_data()
    # h.save_data('/home/mike/research/ac6-microburst-scale-sizes'
    #             '/data/norm/ac6_norm_all.txt')