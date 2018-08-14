# This program will catalogue AC6 microbursts

import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import copy
import sys
import csv
import os

sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
sys.path.insert(0, '/home/mike/research/microburst-detection/burst_parameter')
sys.path.insert(0, '/home/mike/research/microburst-detection/wavelets')
sys.path.insert(0, '/home/mike/research/mission-tools/misc')

import read_ac_data
import burst_parameter
import waveletAnalysis
import locate_consecutive_numbers

class FindMicrobursts(waveletAnalysis.WaveletDetector):
    def __init__(self, sc_id, date, runPipeline=False):
        self.sc_id = sc_id
        self.date = date
        self._loadData() # Read in the data
        if runPipeline: # Run the microburst detector code and save data to CSV file
            self.getMicroburstIdx()
            self.saveData()
        return

    def getMicroburstIdx(self, thresh=5, method='wavelet'):
        """
        Use either obrien or wavelet method to calculate the 
        microbrust indicies _getPeak() method will find the peak
        amplitudes of microbursts
        """
        if method == 'obrien':
            self._getBurstParam()
        else:
            self._getWavelet()

        self._checkMicroburstFlag()
        return

    def saveData(self, fPath=None):
        """
        This function will save microburst data to a CSV file.
        Columns are: dateTime, dos1 amplitude (not baseline subtracted),
        L, MLT, lat, lon, alt, Loss_Cone_Type, 
        """
        keys = ['dateTime', 'dos1rate', 'Lm_OPQ', 'MLT_OPQ', 'lat',
            'lon', 'alt', 'Dist_In_Track', 'Lag_In_Track', 'Dist_Total','Loss_Cone_Type']
        headerl1 = ['Microburst catalogue created on {}'.format(
            datetime.now())]
        headerl2 = copy.copy(keys)
        #headerl2[0] = '# {}'.format(headerl2[0])

        if fPath is None:
            saveDir = os.path.abspath('./../data/z_daily_microburst_catalogues/')
            saveName = 'AC6{}_{}_microbursts.txt'.format(self.sc_id, self.date.date())
            fPath = os.path.join(saveDir, saveName)

        with open(fPath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headerl1)
            writer.writerow(headerl2)

            line = [None]*len(keys)
            for peakInd in self.peakInd:
                for ic, key in enumerate(keys):
                    line[ic] = self.d[key][peakInd]
                writer.writerow(line)
        return

    def _checkMicroburstFlag(self, gTx=True, scTx=True):
        """
        Filter the microburst indicies by the data quality flag.
        Refer to AC-6 data README for more information on the flag
        bitmap.
        
        OPTIONAL ARGS:
        gTx = True : Filter out detections made during ground transmission
        scTx = True : Filter out detections made during crosslink transmission
        """
        flag = np.array(self.d['flag'], dtype=int)
        if gTx:
            # Identify flags without the ground transmission.
            gTxInd = np.bitwise_and(flag, 1)
        else:
            gTxInd = np.zeros_like(flag)
        if scTx:
            # Identify flags without the cross-spacecraft transmission.
            scTxInd = np.bitwise_and(flag, 2)
        else:
            scTxInd = np.zeros_like(flag)
        # Index array of bad flags.
        badBool = np.logical_or(gTxInd, scTxInd) 
        goodBool = np.logical_not(badBool)    
        goodInd = np.where(goodBool)[0]

        sameInd = list(set(self.peakInd) & set(goodInd))
        self.peakInd = np.array(sorted(sameInd), dtype=int)    
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

    def _getBurstParam(self, ch='dos1rate', n=0.1, a=0.3):
        """
        Calculate the burst parameter on the day. This is a CPU intensive task so
        parallize it?
        """
        self.burstParam = burst_parameter.obrien_burst_param(
            self.d[ch], 0.1, N_WIDTH=n, A_WIDTH=a)
        return

    def _getWavelet(self, ch='dos1rate', thresh=0.1):
        """
        This function handles the creation and manipulation of the wavelet
        detector.
        """
        # Feed the counts into the wavelet microburst finder
        validDataIdt = np.where(self.d[ch] != -1E31)[0]
        waveletAnalysis.WaveletDetector.__init__(self, self.d[ch][validDataIdt], 
            self.d['dateTime'][validDataIdt], 0.1, mother='DOG')
        self.waveletTransform() # Get wavelet space
        self.waveletFilter(self.s0, 1) # Do a band pass and significance filter.
        self.degenerateInvWaveletTransform() # Inverse transform filtered data.
        # Indicies where the error-filetered data is greater than thresh
        self.burstIdt = np.where(self.dataFlt > thresh)[0] 
        self._getPeaks(ch, validDataIdt) # Find peaks
        return
        
    def _getPeaks(self, ch, validDataIdt):
        """
        This function will look for periods of consecutive indicies of detections and
        find the peak for each period. The peak index array for the data is self.peakInd.  
        """
        startInd, endInd = locate_consecutive_numbers.locateConsecutiveNumbers(
            self.burstIdt) # Find consecutive numbers to get a max of first
        self.peakInd = np.nan*np.ones(len(startInd), dtype=int)
        # Loop over every microburst detection region (consecutive microburst indicies)
        for i, (st, et) in enumerate(zip(startInd, endInd)):
            # Index nightmare, but works. There may be a better way
            offset = validDataIdt[self.burstIdt[st]]
            self.peakInd[i] = np.argmax(self.d[ch][validDataIdt[self.burstIdt[st:et]]]) + offset
        self.peakInd = self.peakInd.astype(int)
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
        self.ax[0].scatter(self.d['dateTime'][validIdt[self.burstIdt]], self.d['dos1rate'][validIdt[self.burstIdt]], c='b', s=50)
        self.ax[0].scatter(self.d['dateTime'][self.peakInd], self.d['dos1rate'][self.peakInd], c='r', s=25)
        self.ax[1].plot(self.time, self.dataFlt)
        self.ax[0].set_ylim(bottom=0, top=np.max(self.d['dos1rate'][validIdt]))
        self.ax[2].scatter(self.d['dateTime'][self.peakInd], self.d['flag'][self.peakInd], c='r', s=25)
        return

if __name__ == '__main__':
    # Good day for microbursts to test: 2016-10-14
    # Bad day to test: 2015-04-14
    for sc_id in ['A', 'B']:
        date = datetime(2015, 4, 14) 
        obj = TestFindMicrobursts(sc_id, date)
        obj.getMicroburstIdx()
        #obj.saveData()
        obj.plotTimeseries()
        plt.show()
