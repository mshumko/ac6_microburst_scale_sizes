# This program will catalogue AC6 microbursts

import matplotlib.pyplot as plt
import scipy.signal
from datetime import datetime, timedelta
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

    def getMicroburstIdx(self, method='wavelet', **kwargs):
        """
        Use either obrien or wavelet method to calculate the 
        microbrust indicies _getPeak() method will find the peak
        amplitudes of microbursts
        """
        if method == 'obrien':
            self._getBurstParam(**kwargs)
        else:
            self._getWavelet(**kwargs)

        self._checkMicroburstFlag()
        return

    def corrFlag(self, ACwidth=2, CCwidth=1):
        """ 
        This method implements the correlation filters to remove 
        noise.
        """
        tempPeaks = np.array([], dtype=int)
        for p_i in self.peakInd:
            t_i = self.d['dateTime'][p_i]
            autoCorrRange = [t_i-timedelta(seconds=ACwidth/2), 
                            t_i+timedelta(seconds=ACwidth/2)] 
            crossCorrRange = [t_i-timedelta(seconds=CCwidth/2), 
                            t_i+timedelta(seconds=CCwidth/2)]
            validIdA = np.where(
                        (self.d['dateTime'] > autoCorrRange[0]) & 
                        (self.d['dateTime'] < autoCorrRange[1])
                        )[0]
            # Autocorrelate dos1
            ac1, lags1, max1, peakInd1 = self._autoCorrCounts(
                                    self.d['dateTime'], 
                                    self.d['dos1rate'], 
                                    autoCorrRange)
            # Autocorrelate dos2 (because dos1 and dos2 respond 
            # similarly to noise)
            ac2, lags2, max2, peakInd2 = self._autoCorrCounts(
                                    self.d['dateTime'], 
                                    self.d['dos2rate'], 
                                    autoCorrRange)
            # Cross correlate dos1 and dos2
            dos12corr, dos12lags = self._crossCorrCounts(
                                    self.d['dateTime'], 
                                    self.d['dateTime'], 
                                    self.d['dos1rate'], 
                                    self.d['dos2rate'], 
                                    crossCorrRange)                            
            # First check that dos1 and dos12 are correlated, then
            # then if max(dos2) > 0.5*max(dos1) then it is noise. 
            if not ( (len(np.where(peakInd1 == 4)[0]) == 1 or 
                    len(np.where(peakInd1 == 2)[0]) == 1 or 
                    len(np.where(peakInd2 == 4)[0]) == 1 or 
                    len(np.where(peakInd2 == 2)[0]) == 1) and
                    (max(self.d['dos2rate'][validIdA]) > 1000 or
                    max(dos12corr) >= 0.9) ):
                tempPeaks = np.append(tempPeaks, p_i)    
        self.peakInd = tempPeaks
        return


    def saveData(self, fPath=None):
        """
        This function will save microburst data to a CSV file.
        Columns are: dateTime, dos1 amplitude (not baseline subtracted),
        L, MLT, lat, lon, alt, Loss_Cone_Type, 
        """
        keys = ['dateTime', 'dos1rate', 'Lm_OPQ', 'MLT_OPQ', 'lat',
            'lon', 'alt', 'Dist_In_Track', 'Lag_In_Track',
            'Dist_Total','Loss_Cone_Type', 'flag']
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

    def _autoCorrCounts(self, times, counts, tRange, norm=True):
        """
        This method calculates the autocorrelation of the counts array in
        the time range specified by tRange. times array is used to identify 
        which counts to autocorrelate.
        
        The ax argument specified the subplot on which to plot the 
        autocorrelation on. 
        """
        validIdt = np.where((counts != -1E31) & 
                            (times > tRange[0]) & 
                            (times < tRange[1]))[0]
        x = counts[validIdt] - counts[validIdt].mean()
        # mode=same means that some edge effects will be observed. Should be ok.                    
        ac = np.correlate(x, x, mode='same')
        ac = ac[ac.size//2:]
        # Lags are normalized to seconds  
        lags = np.arange(0, ac.size)/10 
        
        if norm:
            ac /= len(x)*np.var(x)
            
        # Identify peaks
        peakInd, _ = scipy.signal.find_peaks(ac, prominence=0.1)
        return ac, lags, max(counts[validIdt]), peakInd

    def _crossCorrCounts(self, timesA, timesB, countsA, countsB, tRange, norm=True):
        """
        This method calculates the crosscorrelation of the counts array in
        the time range specified by tRange. times array is used to identify 
        which counts to autocorrelate. This method also calculates the cross
        correlation of the two signals to determine if the micorbursts 
        observed by both AC6 units are coincident.
        
        The ax argument specified the subplot on which to plot the 
        autocorrelation on. 
        """
        validIdtA = np.where((countsA != -1E31) & 
                            (timesA > tRange[0]) & 
                            (timesA < tRange[1]))[0]
        x = countsA[validIdtA] - countsA[validIdtA].mean()

        validIdtB = np.where((countsB != -1E31) & 
                            (timesB > tRange[0]) & 
                            (timesB < tRange[1]))[0]
        y = countsB[validIdtB] - countsB[validIdtB].mean()

        cc = np.correlate(x, y, mode='same')
        # Lags are normalized to seconds  
        lags = np.arange(-cc.size/2, cc.size/2)/10 
        
        if norm:
            cc /= np.sqrt(len(x)*np.var(x)*len(y)*np.var(y))
            
        # Identify peaks
        #peakInd, _ = scipy.signal.find_peaks(cc)
        return cc, lags

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

    def _getBurstParam(self, ch='dos1rate', n=0.1, a=0.5, thresh=5):
        """
        Calculate the burst parameter on the day. This is a CPU intensive task so
        parallize it?
        """
        validDataIdt = np.where(self.d[ch] != -1E31)[0]
        self.burstParam = burst_parameter.obrien_burst_param(
            self.d[ch][validDataIdt], 0.1, N_WIDTH=n, A_WIDTH=a)
        self.burstIdt = np.where(self.burstParam > thresh)[0] 
        self._getPeaks(ch, validDataIdt) # Find peaks
        return

    def _getWavelet(self, ch='dos1rate', thresh=0.1, maxWidth=1, SIGNIF_LEVEL=0.25):
        """
        This function handles the creation and manipulation of the wavelet
        detector.
        """
        # Feed the counts into the wavelet microburst finder
        validDataIdt = np.where(self.d[ch] != -1E31)[0]
        waveletAnalysis.WaveletDetector.__init__(self, self.d[ch][validDataIdt], 
            self.d['dateTime'][validDataIdt], 0.1, mother='DOG', siglvl=0.95)
        self.waveletTransform() # Get wavelet space
        self.waveletFilter(self.s0, maxWidth, SIGNIF_LEVEL=SIGNIF_LEVEL) # Do a band pass and significance filter.
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
            if st == et: 
                # If the same index 
                et += 1
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
        self.ax[0].plot(self.d['dateTime'][validIdt], 
                        self.d['dos1rate'][validIdt])
        self.ax[0].fill_between(self.d['dateTime'][validIdt], 
            self.d['dos1rate'][validIdt]-np.sqrt(self.d['dos1rate'][validIdt]),
            self.d['dos1rate'][validIdt]+np.sqrt(self.d['dos1rate'][validIdt]),
            color='r', alpha=0.5)
        self.ax[0].scatter(self.d['dateTime'][validIdt[self.burstIdt]], self.d['dos1rate'][validIdt[self.burstIdt]], c='b', s=50)
        self.ax[0].scatter(self.d['dateTime'][self.peakInd], self.d['dos1rate'][self.peakInd], c='r', s=25)
        if hasattr(self, 'time'):
            self.ax[1].plot(self.d['dateTime'][validIdt], self.dataFlt)
        else:
            self.ax[1].plot(self.d['dateTime'][validIdt], self.burstParam)
            self.ax[1].set_ylim(-10, 10)
        self.ax[0].set_ylim(bottom=0, top=np.max(self.d['dos1rate'][validIdt]))
        self.ax[2].scatter(self.d['dateTime'][self.peakInd], self.d['flag'][self.peakInd], c='r', s=25)
        return

if __name__ == '__main__':
    # Good day for microbursts to test: 2016-10-14
    # Bad day to test: 2015-04-14
    for sc_id in ['A', 'B']:
        date = datetime(2016, 10, 14) 
        obj = TestFindMicrobursts(sc_id, date)
        obj.getMicroburstIdx(method='obrien')
        #obj.saveData()
        obj.plotTimeseries()
        plt.show()
