import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import sys

sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

class SpaceTime:
    def __init__(self, date):
        """
        This class plots the AC-6 time series aligned in time 
        (by time stamps), and in space by accounting for the 
        in-track lag.
        """
        self.date = date
        self._load_data(date)

        return

    def plot(self, ax=None, tRange=None):
        """
        This method plots the AC-6 data in time aligned, 
        and track-aligned forms.
        """
        if ax is None: 
            _, self.ax = plt.subplots(2, sharex=True, figsize=(12, 8))
        else:
            self.ax = ax

        if tRange is not None:
            # Find filtered indicies
            idta = np.where((self.aData['dateTime'] > tRange[0]) & 
                            (self.aData['dateTime'] < tRange[1]) &
                            (self.aData['dos1rate'] != -1E31))[0]
            idtb = np.where((self.bData['dateTime'] > tRange[0]) & 
                            (self.bData['dateTime'] < tRange[1]) &
                            (self.bData['dos1rate'] != -1E31))[0]
            # idtbShifted = np.where((self.bData['dateTime'] > tRange[0] + 
            #                 np.mean(self.bData['Lag_In_Track'])) & 
            #                 (self.bData['dateTime'] < tRange[1] + 
            #                  np.mean(self.bData['Lag_In_Track']))) &
            #                 (self.bData['dos1rate'] != -1E31))[0]
            # filter data
            filterKeys = ['dateTime', 'dos1rate', 'dos2rate', 
                          'dos3rate', 'Lm_OPQ', 'MLT_OPQ']
            tempDictA = {key:self.aData[key][idta] for key in filterKeys}
            tempDictB = {key:self.bData[key][idtb] for key in filterKeys}
        else:
            tempDictA = self.aData
            tempDictB = self.bData

        # 
        self._time_plot(self.ax[0], tempDictA['dateTime'], 
                        tempDictB['dateTime'], 
                        tempDictA['dos1rate'], 
                        tempDictB['dos1rate'])
        self._space_plot(self.ax[1], tRange)
        
        return

    def _time_plot(self, ax, timeA, timeB, dosA, dosB):
        """
        This methods plots the time-aligned time series.
        """
        ax.plot(timeA, dosA, 'r', label='AC6-A')
        ax.plot(timeB, dosB, 'b', label='AC6-B')
        ax.legend(loc=1)
        ax.set(ylabel='dos rate [counts/s]', 
              yscale='log',
              title=('AC-6 Space-Time plot | '
                    '{}'.format(self.date.date())))
        return

    def _space_plot(self, ax, tRange):
        """
        This methods plots the space-aligned time series.
        """
        # Shift AC6B times
        # Find the in-track lag at this time first
        idt = np.where(self.bData['dateTime'] > tRange[0])[0] 
        lag = self.bData['Lag_In_Track'][idt[0]]
        bTimes = np.array([t - timedelta(seconds=lag) for t in self.bData['dateTime']])
        
        # Find the valid indicies to plot.
        idta = np.where((self.aData['dateTime'] > tRange[0]) & 
                            (self.aData['dateTime'] < tRange[1]) &
                            (self.aData['dos1rate'] != -1E31))[0]
        idtb = np.where((bTimes > tRange[0]) & 
                        (bTimes < tRange[1]) &
                        (self.bData['dos1rate'] != -1E31))[0]
        ax.plot(self.aData['dateTime'][idta], self.aData['dos1rate'][idta], 'r')
        ax.plot(bTimes[idtb], self.bData['dos1rate'][idtb], 'b')
        ax.set(ylabel='dos rate [counts/s]', 
              yscale='log', xlabel='UTC')
        ax.text(0.05, 1, 'Dist_In_Track = {} km'
                '\n Lag_In_Track = {} s'.format(
                round(self.bData['Dist_In_Track'][idt[0]]), 
                round(lag, 1)), va='top', ha='left', transform=ax.transAxes)
        return


    def _load_data(self, date):
        """
        This method loads in both of the AC6 datasets.
        """        
        self.aData = read_ac_data.read_ac_data_wrapper('A', date,
            dType='10Hz', plot=False)
        self.bData = read_ac_data.read_ac_data_wrapper('B', date,
            dType='10Hz', plot=False)
        return

if __name__ == '__main__':
    tRange = [datetime(2016, 10, 14, 4, 26), datetime(2016, 10, 14, 4, 29)]
    #tRange = [datetime(2017, 1, 11, 1, 52, 27), datetime(2017, 1, 11, 1, 54)]
    st = SpaceTime(tRange[0])
    st.plot(tRange=tRange)
    plt.show()