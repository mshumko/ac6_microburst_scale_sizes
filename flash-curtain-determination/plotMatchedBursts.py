import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates
import sys
import os
import scipy.signal
from datetime import datetime, timedelta

sys.path.append('/home/mike/research/ac6-microburst-scale-sizes/stats/')
sys.path.append('/home/mike/research/mission-tools/ac6/')
import read_ac_data
import scale_sizes

class PlotMicroburstMatches(scale_sizes.ScaleSize):
    def __init__(self, cPath, burstType='flashes'):
        """
        This class plots the matched microbursts for AC-6, flashes and 
        curtains.
        """
        scale_sizes.ScaleSize.__init__(self, cPath) # Read in catalogue data.

        # Shift times if want to plot curtains
        if 'curtain' in burstType:
            # shift times

            # CHECK IF THIS SHIFT NEEDS TO BE DONE HERE OR WHEN PLOTTING!
            self.dataDict['dateTimeB'] = np.array([t + timedelta(seconds=dt) 
                for (t, dt) in zip(self.dataDict['dateTimeB'], 
                self.dataDict['Lag_In_Track'])])
            self.burstType = 'curtain'
        elif 'flash' in burstType:
            self.burstType = 'flash'
        else:
            raise NameError("burstType kwarg must be either 'flash' or 'curtain'.")
        return

    def plotMicroburst(self, thresh=2):
        """ 
        This function will plot the microburst flash or curtain detections.
        """
        self._findDateBounds() # Get time bounds for unique dates

        fig, self.ax = plt.subplots()
        self.bx = self.ax.twinx()

        # Loop over dates.
        for d in self.dateArr:
            # Open file
            d[0] = datetime.combine(d[0], datetime.min.time())
            try: # Remove the try except code when I have access to all of the data.
                self.dataA = read_ac_data.read_ac_data_wrapper('A', d[0])
                self.dataB = read_ac_data.read_ac_data_wrapper('B', d[0])
            except AssertionError as err:
                if 'None or > 1 AC6 files found in' in str(err):
                    continue
                else:
                    raise

            # Loop over daily detections.
            zippedCenterTimes = zip(self.dataDict['dateTimeA'][d[1]:d[2]], 
                self.dataDict['dateTimeB'][d[1]:d[2]])
            for (tA, tB) in zippedCenterTimes: 
                # Plot dos1rate
                pltRange = [tA-timedelta(seconds=thresh), 
                            tA+timedelta(seconds=thresh)]
                flag = self.plotTimeRange(self.dataA['dateTime'], 
                    self.dataA['dos1rate'],
                    self.dataB['dateTime'], self.dataB['dos1rate'], pltRange)

                # Save figure
                saveDir = ('/home/mike/research/ac6-microburst-scale-sizes/'
                    'data/plots/validation/{}/'.format(datetime.now().date()))
                if not os.path.exists(saveDir): # Create save dir if it does not exist.
                    os.makedirs(saveDir)

                saveDate = pltRange[0].replace(microsecond=0).isoformat().replace(':', '')
                saveName = '{}_validation_{}.png'.format(self.burstType, saveDate)
                # Format times and plots.
                tfmt = matplotlib.dates.DateFormatter('%H:%M:%S')
                self.ax.xaxis.set_major_formatter(tfmt)
                plt.tight_layout()
                if flag == 1:
                    print('Saving figure {}'.format(saveName))
                    plt.savefig(os.path.join(saveDir, saveName))     
                self.ax.clear()
                self.bx.clear()
        return

    def plotTimeRange(self, timeA, dosA, timeB, dosB, tRange):
        """
        This function will plot the narrow microburst time range from
        both spacecraft, and annotate with L, MLT, Lat, Lon. If curtain,
        add shift in time and space.
        """
        validIdA = np.where((dosA != -1E31) & (timeA > tRange[0]) & (timeA < tRange[1]))[0]
        validIdB = np.where((dosB != -1E31) & (timeB > tRange[0]) & (timeB < tRange[1]))[0]

        # Need to pass in the time array as well. Maybe do the time shift here?
        self.ax.plot(timeA[validIdA], dosA[validIdA], 'r', label='AC-6 A')
        self.ax.plot(timeB[validIdB], dosB[validIdB], 'b', label='AC-6 B')
        self.bx.plot(timeA[validIdA], self.dataA['Alpha'][validIdA], '--r')
        self.bx.plot(timeB[validIdB], self.dataB['Alpha'][validIdB], '--b')
        self.ax.set(ylabel='Dos1 [counts/s]', xlabel='UTC',
            title='AC-6 {} validation | {}'.format(
            self.burstType, tRange[0].date()))
        self.bx.set(xlabel='UTC', ylabel='Alpha [Deg] (dashed)')
        self.ax.legend(loc=1)

        # Do the cross correlation
        cc = self._calcCrossCorr(dosA, dosB, validIdA, validIdB)
        
        if cc < 0.9: # Simple cross correlation test
            return 0
        
        # Calc power spectrum
#        f, ps = self._calcFrequSpec(dosA[validIdA])
#        self.ax.plot(f, ps, 'r')
#        f, ps = self._calcFrequSpec(dosB[validIdB])
#        self.ax[.plot(f, ps, 'b')
#        self.ax[1].set_xlim(left=0)
        meanFlag = (np.mean(self.dataA['flag'][validIdA]) + 
            np.mean(self.dataB['flag'][validIdB]))/2
        textStr = ('L={} MLT={}\nlat={} lon={}\ndist={} LCT={}\nCrossCorr={} flag={}'.format(
            round(np.mean(self.dataA['Lm_OPQ'][validIdA])),
            round(np.mean(self.dataA['MLT_OPQ'][validIdA])),
            round(np.mean(self.dataA['lat'][validIdA])),
            round(np.mean(self.dataA['lon'][validIdA])),
            round(np.mean(self.dataA['Dist_In_Track'][validIdA])),
            round(np.mean(self.dataA['Loss_Cone_Type'][validIdA])),
            round(cc, 1), round(meanFlag)))
        self.ax.text(0.05, 0.95, textStr, transform=self.ax.transAxes, 
            va='top')
        return 1
        
    def _calcCrossCorr(self, dosA, dosB, validIdtA, validIdtB, width=5):
        """
        This function calculates a cross correlation of the data with a
        data width in data points. 
        """
        corrIndA = validIdtA[len(validIdtA)//2-width : len(validIdtA)//2+width]
        corrIndB = validIdtB[len(validIdtB)//2-width : len(validIdtB)//2+width]
        cc = scipy.signal.correlate(dosA[corrIndA] - np.mean(dosA[corrIndA]), 
            dosB[corrIndB] - np.mean(dosB[corrIndB]), mode='valid')[0]
        cc /= np.std(dosA[corrIndA])*np.std(dosB[corrIndB])*(len(corrIndA)+len(corrIndB))/2
        return cc
        
    def _calcFrequSpec(self, dosRate, cadence=0.1):
        """
        This function calculates the frequency spectrum on the dosRate data
        using an FFT with a 10 Hz data by default. Inspired by 
        https://stackoverflow.com/questions/15382076/plotting-power-spectrum-in-python
        """
        dosRate -= np.mean(dosRate) # Remove DC offset.
        ps = np.abs(np.fft.fft(dosRate))**2
        f = np.fft.fftfreq(len(dosRate), cadence)
        idx = np.argsort(f) # Not sure why these need to be sorted, but cant hurt.
        return f[idx], ps[idx]


    def _findDateBounds(self):
        """ 
        This function will calculate the unique dates in the catalogue file
        to speed up the plotting reutine.
        """
        dates = np.array([t.date() for t in self.dataDict['dateTimeA']])
        uniqueDates = np.array(sorted(set(dates)))
        #print(uniqueDates)

        # Make an array with columns with the unique date, startInd, endInd.
        self.dateArr = np.nan*np.ones((len(uniqueDates), 3), dtype=object)
        # Fill in the self.dateArr with the start/end indicies that correspond
        # to the catalogue.
        for (i, d) in enumerate(uniqueDates):
            #print(dates, d)
            self.dateArr[i, 0] = d
            dateInd = np.where(dates == d)[0]
            #print(d, dateInd)
            self.dateArr[i, 1] = dateInd[0]
            self.dateArr[i, 2] = dateInd[-1]
        return

if __name__ == '__main__':
    cPath = ('/home/mike/research/ac6-microburst-scale-sizes/data/'
        'flash_catalogues/flashes_catalogue_v1.txt')
    pltObj = PlotMicroburstMatches(cPath)
    pltObj.plotMicroburst()
