import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
import datetime

import sys
sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

class OccuranceRate:
    def __init__(self, sc_id, date):
        """
        For a given spacecraft on a given day, this class calculates
        the microburst occurance rate (with an assumption about width.)
        """
        self.sc_id = sc_id
        self.date = date
        self._load_data()
        return

    def radBeltIntervals(self, lbound=4, ubound=8):
        """
        This method finds all of the radiation belt passes in a day 
        and separates it into start/end indicies.
        """
        L = np.abs(self.data['Lm_OPQ'])
        idL = np.where((L > lbound) & (L < ubound))[0]

        conv = np.convolve([1, -1], idL, mode = 'valid') - 1
        consecutiveFlag = np.where(conv != 0)[0] + 1
        startInd = np.insert(consecutiveFlag, 0, 0)
        endInd = np.insert(consecutiveFlag, 
                len(consecutiveFlag), len(idL)-1)
        
        self.intervals = np.zeros((len(startInd), 2), dtype=int)
        for i, (i_s, i_e) in enumerate(zip(startInd, endInd)):
            self.intervals[i, :] = [idL[i_s], idL[i_e-1]]
        return 

    def _load_data(self):
        """ Loads AC-6 10 Hz data """
        self.data = read_ac_data.read_ac_data_wrapper(self.sc_id, self.date)
        return

if __name__ == '__main__':
    o = OccuranceRate('A', datetime.datetime(2016, 10, 14))
    o.radBeltIntervals()

    plt.plot(o.data['dateTime'], o.data['Lm_OPQ'], 'b')

    # for (sInd, eInd) in zip(o.startInd, o.endInd):
    #     plt.scatter(np.arange(sInd, eInd), 
    #                 o.data['Lm_OPQ'][sInd:eInd], c='r') 

    for (sInd, eInd) in o.intervals:
        print(sInd, eInd)
        plt.scatter(o.data['dateTime'][sInd:eInd], 
                    o.data['Lm_OPQ'][sInd:eInd], c='r') 


    plt.ylim(0, 12)
    plt.show()
