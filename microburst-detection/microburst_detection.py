# This program will catalogue AC6 microbursts

import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import sys
import os

sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
sys.path.insert(0, '/home/mike/research/microburst-detection/burst_parameter')

import read_ac_data
import burst_parameter

class FindMicrobursts:
    def __init__(self, sc_id, date):
        self.sc_id = sc_id
        self.date = date
        self._loadData() # Read in the data
        self.getMicroburstIdx() # Calculate microburst indicies
        return

    def getMicroburstIdx(self, thresh=10, method='obrien'):
        """
        Use either obrien or wavelet method to calculate the 
        microbrust indicies
        """
        if method == 'obrien':
            self._getBurstParam()
        else:
            raise NotImplemented('Wavelet method not implemented yet!')

        self.burstIdt = np.where(self.burstParam > thresh)[0]
        #self._checkMicroburstFlag()
        return

    def _checkMicroburstFlag(self):
        """
        Filter the microburst indicies by the data quality flag.
        """
        self.validFlagIdt = (self.d['flag'] == 0)
        sameIdt = list(set(self.burstIdt) & set(self.validFlagIdt))
        self.burstIdt = np.array(sorted(sameIdt))        
        return

    def _loadData(self):
        """
        Load the AC-6 data.

        10 Hz data keys:
            'alt', 'lat', 'lon', 'dos1l', 'dos1m', 'dos1rate', 'dos2l', 'dos2m',
            'dos2rate', 'dos3l', 'dos3m', 'dos3rate', 'flag', 'Subcom', 'Lm_OPQ', 
            'Bmag_OPQ', 'MLT_OPQ', 'InvLat_OPQ', 'Loss_Cone_Type', 'K_Z', 'Lstar_Z',
            'hmin_Z', 'Alpha', 'Beta', 'Dist_In_Track', 'Lag_In_Track', 
            'Dist_Cross_Track_Horiz', 'Dist_Cross_Track_Vert', 'Dist_Total'
        """
        self.d = read_ac_data.read_ac_data_wrapper(self.sc_id, self.date,
            dType='10Hz', plot=False)
        return

    def _getBurstParam(self, ch='dos1rate', n=0.1, a=0.5):
        """
        Calculate the burst parameter on the day. This is a CPU intensive task so
        parallize it?
        """
        self.burstParam = burst_parameter.obrien_burst_param(
            self.d[ch], 0.1, N_WIDTH=n, A_WIDTH=a)
        return


class TestFindMicrobursts(FindMicrobursts):
    def __init__(self, sc_id, date):
        self.sc_id = sc_id
        self.date = date
        FindMicrobursts.__init__(self, self.sc_id, self.date)  

        # Create empty plots
        self.fig, self.ax = plt.subplots(3, sharex=True)
        return

    def plotTimeseries(self):
        validIdt = np.where(self.d['dos1rate'] != -1E31)[0]
        self.ax[0].plot(self.d['dateTime'][validIdt], self.d['dos1rate'][validIdt])
        self.ax[1].plot(self.d['dateTime'], self.d['flag'])
        self.ax[2].plot(self.d['dateTime'][self.burstIdt], self.d['dos1rate'][self.burstIdt])
        return
    


if __name__ == '__main__':
    sc_id = 'A'
    date = datetime(2016, 9, 30)
    obj = TestFindMicrobursts(sc_id, date)
    obj.plotTimeseries()
    plt.show()