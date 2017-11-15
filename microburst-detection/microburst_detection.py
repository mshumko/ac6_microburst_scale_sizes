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
        self._loadData()
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
            dType='10Hz', plot=True)
        return


class TestFindMicrobursts(FindMicrobursts):
    def __init__(self, sc_id, date):
        self.sc_id = sc_id
        self.date = date
        FindMicrobursts.__init__(self, self.sc_id, self.date)  
        return


if __name__ == '__main__':
    sc_id = 'A'
    date = datetime(2016, 9, 30)
    obj = TestFindMicrobursts(sc_id, date)
    plt.plot()