# This function calculates the normalization parameters for the 
# micorburst scale size distributions.

from datetime import datetime, timedelta
from matplotlib.dates import date2num 
import csv
import numpy as np
import sys
import os

# Import personal libraries
sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

class Hist1D:
    def __init__(self, d=None, startDate=datetime(2014, 1, 1),
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
        self.filterDict = filterDict
        self.flag = flag
        return

    def loop_data(self):
        """

        """
        for day in self.dates:
            self.load_day_data(day)
            if (self.ac6dataA is None) or (self.ac6dataB is None):
                continue # If one (or both) of the data files is empty
            ind = self.filterData()
            # If using Hist2D, the Hist1D's method will be overwritten.
            self.hist_data(ind) 

        return

    def load_day_data(self, date):
        """
        This generator function will load in one day of data at a 
        time, and return control to user.
        """
        try:   
            print('Loading data from {}'.format(date))
            self.ac6dataA = read_ac_data.read_ac_data_wrapper(
                'A', date, dType='10Hz', plot=False)
            self.ac6dataB = read_ac_data.read_ac_data_wrapper(
                'B', date, dType='10Hz', plot=False)
        except AssertionError as err:
            if ( ('None or > 1 AC6 files found' in str(err)) or
                ('File is empty!' in str(err)) ):
                self.ac6dataA = None
                self.ac6dataB = None
            else:
                raise
        finally:
            return self.ac6dataA, self.ac6dataB

    def filterData(self, verbose=True):
        """
        This function filters the AC-6 data by common times, data flag value,
        and filterDict dictionary.
        """
        if verbose:
            start_time = datetime.now()
            print('Filtering data at {}'.format(datetime.now()))
        ### Find common times of the two data sets ###
        tA = date2num(self.ac6dataA['dateTime'])
        tB = date2num(self.ac6dataB['dateTime'])
        # np.in1d returns a boolean array that correspond to indicies in tB
        # that are also in tA. np.where will convert this mask array into
        # an index array
        ind = np.where(np.in1d(tB, tA, assume_unique=True))[0]

        ### Data quality flag filter ###
        if self.flag: # First filter by common times and flag
            indf = np.where(self.ac6dataB['flag'] == 0)[0]
            ind = np.intersect1d(ind, indf) 

        ### filerDict filter ###
        for key, value in self.filterDict.items():
            idx = np.where((self.ac6dataB[key] > np.min(value)) & 
                            (self.ac6dataB[key] < np.max(value)))[0]
            ind = np.intersect1d(ind, idx)
        if verbose:
            print('Data filted in {} s'.format(
                (datetime.now()-start_time).total_seconds()))
        return ind

    def hist_data(self, ind):
        """
        This method will histrogram the total distance data.
        """
        H, _ = np.histogram(self.ac6dataB['Dist_Total'][ind], bins=self.d)
        self.count += H/10
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
    def __init__(self, histKeyX, histKeyY, bins=None, startDate=datetime(2014, 1, 1),
                 endDate=datetime.now(), filterDict={}, flag=True):
        """
        This class calculates the 2D histograms as a function of distance
        for the various filter parameters. 
        """
        Hist1D.__init__(self, d=None, startDate=startDate, 
                        endDate=endDate, filterDict=filterDict, flag=flag)
        self.histKeyX = histKeyX
        self.histKeyY = histKeyY

        if bins is None: # Determine what the bins should be.
            if 'lm' in histKeyX.lower(): # L-MLT distribution
                self.bins=np.array([np.arange(2, 10), np.arange(0, 24)])
            elif 'lm' in histKeyY.lower():
                self.bins=np.array([np.arange(0, 24), np.arange(2, 10)])

            elif 'lat' in histKeyX.lower(): # lat-lon distribution
                self.bins=np.array([np.arange(-90, 90, 5), np.arange(-180, 180, 5)])
            elif 'lat' in histKeyY.lower():
                self.bins=np.array([np.arange(-90, 90, 5), np.arange(-180, 180, 5)])
            else:
                raise ValueError('Could not determine how to make bins!'
                                ' Plase supply them')
        else:
            self.bins = bins
        
        self.count = np.zeros((len(self.bins[0])-1, len(self.bins[1])-1))
        return

    def hist_data(self, ind):
        """
        This histogram method overwrites the Hist1D's hist_data() method
        """
        H, xedges, yedges = np.histogram2d(self.ac6dataB[self.histKeyX][ind],
                                self.ac6dataB[self.histKeyY][ind], bins=self.bins)

        self.count += H/10
        return

    def save_data(self, fPathBin, fPathNorm):
        """
        This method saves the histrogram bins and normalization coefficients.
        """
        XX, YY = np.meshgrid(self.bins[0], self.bins[1])

        with open(fPathBin, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow([self.histKeyX, self.histKeyY])
            w.writerows(self.bins)

        with open(fPathNorm, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow([self.histKeyX, self.histKeyY])
            w.writerows(self.count)
        return


if __name__ == '__main__':
    ### SCRIPT TO MAKE "Dst_Total" NORMALIZATION ###
    # ss=Hist1D()
    # st = datetime.now()
    # ss.loop_data()
    # sDir = '/home/mike/research/ac6-microburst-scale-sizes/data/norm/'
    # ss.save_data(os.path.join(sDir, 'ac6_norm_all_test.csv'))
    # print('Norm.py ran in :{} s'.format((datetime.now()-st).total_seconds()))

    ### SCRIPT TO MAKE L-dependent "Dst_Total" NORMALIZATION ###
    # st = datetime.now()
    # L = [3, 4, 5, 6, 7]
    # for (lL, uL) in zip(L[:-1], L[1:]):
    #     ss=Hist1D(filterDict={'Lm_OPQ':[lL, uL]})
    #     ss.loop_data()
    #     sDir = '/home/mike/research/ac6-microburst-scale-sizes/data/norm/'
    #     ss.save_data(os.path.join(sDir, 'ac6_norm_{}_L_{}.csv'.format(lL, uL)))
    # print('Norm.py ran in :{} s'.format((datetime.now()-st).total_seconds()))

    ### SCRIPT TO MAKE L-MLT NORMALIATION ###
    ss = Hist2D('Lm_OPQ', 'lon', bins=[np.arange(2, 10), np.arange(-180, 181, 5)])
    ss.loop_data()
    sDir = '/home/mike/research/ac6-microburst-scale-sizes/data/norm/'
    ss.save_data(os.path.join(sDir, 'ac6_L_lon_bins.csv'), 
                 os.path.join(sDir, 'ac6_L_lon_norm.csv'))

    ss2 = Hist2D('Lm_OPQ', 'MLT_OPQ', bins=[np.arange(2, 10), np.arange(0, 25)])
    ss2.loop_data()
    ss2.save_data(os.path.join(sDir, 'ac6_L_MLT_bins.csv'), 
                  os.path.join(sDir, 'ac6_L_MLT_norm.csv'))