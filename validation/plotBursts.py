import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates
import sys
import os
import itertools
import dateutil.parser
import csv
import scipy.signal
from datetime import datetime, timedelta

#sys.path.append('/home/mike/research/ac6-microburst-scale-sizes/stats/')
sys.path.append('/home/mike/research/mission-tools/ac6/')
import read_ac_data
#import leo_scale_sizes

class ValidateDetections:
    def __init__(self, sc_id, cPathPrimary, cPathSecondary, saveDir=None,
                 uniqueFolders=True):
        """
        This class contains code to validate the microburst detector algorithm
        by plotting the detections made from each spacecraft side by side. No
        attempt is made to match coincident events. 

        For a user-suplied catalog file (path specified by cPath), the 
        catalog file is first read. Then 

        uniqueFolders flag determines if the true/false plots are to be 
        sorted into different folders to easier inspection.
        """
        # Load catalog
        self.cDataPrimary = self._load_catalog_(cPathPrimary)
        self.cDataSecondary = self._load_catalog_(cPathSecondary)
        self.primary_sc = sc_id
        self.uniqueFolders = uniqueFolders
        self.sep_range = [-2E31, 1000]
        # Set save directory
        if saveDir is not None:
            self.saveDir = saveDir
        elif saveDir is None and not uniqueFolders:
            self.saveDir = ('/home/mike/research/ac6-microburst-scale-sizes/'
                    'plots/validation/{}/'.format(datetime.now().date()))
            if not os.path.exists(self.saveDir):
                print('Creating plot directory at:', self.saveDir)
                os.makedirs(self.saveDir)
        else:
            self.saveDir = ('/home/mike/research/ac6-microburst-scale-sizes/'
                    'plots/validation/{}/'.format(datetime.now().date()))
            if not os.path.exists(os.path.join(self.saveDir, 'true')):
                os.makedirs(os.path.join(self.saveDir, 'true'))
            if not os.path.exists(os.path.join(self.saveDir, 'false')):
                os.makedirs(os.path.join(self.saveDir, 'false'))
                print('Created plot directories:', 
                    os.path.join(self.saveDir, 'true'), 'and', 
                    os.path.join(self.saveDir, 'false'))
        return
        
    def plotLoop(self, pltWidth=4, autoCorrWidth=2, crossCorrWidth=1):
        """
        This method loops over the detections in the catalog file and plots 
        the dos1 data.

        The autoCorrWidth and crossCorrWidth kwargs set the correlation window
        sizes for the experimental filtering algirthms. maxDist kwarg sets the
        maximum distance that will be plotted (usefull if you want to plot events
        when the spacecraft were near each other).
        """
        #currentDate = datetime.date().min
        self._findDateBounds() # Get index bounds for each date.
        
        fig, self.ax = plt.subplots(4, figsize=(8, 10))
        
        # Loop over dates.
        for d in self.dateArr:
        
            # If that day's mean separation was outside the self.sep_range 
            # range, move on
            medianDist = np.median(self.cDataPrimary['Dist_Total'][d[1]:d[2]])
            if medianDist < self.sep_range[0] or medianDist > self.sep_range[1]: 
                continue
        
            print('Making validation plots from:', d[0])
            # Open file
            d[0] = datetime.combine(d[0], datetime.min.time())
            try:
                self.dataA = read_ac_data.read_ac_data_wrapper('A', d[0])
                self.dataB = read_ac_data.read_ac_data_wrapper('B', d[0])
            except AssertionError as err:
                if (('None or > 1 AC6 files found in' in str(err)) or 
                    'File is empty!' in str(err)):
                    continue
                else:
                    raise
                
            # The following block of code finds the detections made on day d[0].
            primaryDetTimes = self.cDataPrimary['dateTime'][d[1]:d[2]]
            secondaryDates = np.array([t.date() for t in self.cDataSecondary['dateTime']])
            secondaryIdx = np.where(secondaryDates == d[0].date())[0]
            secondaryDetTimes = self.cDataSecondary['dateTime'][secondaryIdx]

            # Loop over detections on d'th day.
            for i, tA in enumerate(self.cDataPrimary['dateTime'][d[1]:d[2]]):
                # This if statement filters out events detected at 
                # separations > 100 km.
                #if self.cDataPrimary['Dist_Total'][d[1]+i] > maxDist:
                #    continue
                
                pltRange = [tA-timedelta(seconds=pltWidth), 
                            tA+timedelta(seconds=pltWidth)] 
                # Do not plot anything other then good data (flag = 0)
                validIdA = np.where((self.dataA['dateTime'] > pltRange[0]) & 
                                    (self.dataA['dateTime'] < pltRange[1]))[0]
                #if np.any(self.dataA['flag'][validIdA]):
                #    continue
                # Plot dos1rate, centered on the detection.
                flag = self.plotCounts(self.ax[0], self.dataA['dateTime'], 
                    self.dataB['dateTime'], 'dos1rate', pltRange, 
                    primaryDetTimes, secondaryDetTimes)
                if flag == -1:
                    continue

                # Plot dos2rate for the same time period
                flag = self.plotCounts(self.ax[1], self.dataA['dateTime'], 
                    self.dataB['dateTime'], 'dos2rate', pltRange, [], [])

                # Plot the auto-correlation and cross-correlations.
                autoCorrRange = [tA-timedelta(seconds=autoCorrWidth/2), 
                            tA+timedelta(seconds=autoCorrWidth/2)] 
                crossCorrRange = [tA-timedelta(seconds=crossCorrWidth/2), 
                            tA+timedelta(seconds=crossCorrWidth/2)]     
                ac1, lags1, max1, peakInd1 = self.autoCorrCounts(self.dataA['dateTime'], 
                                    self.dataA['dos1rate'], 
                                    autoCorrRange, ax=self.ax[-2])
                ac2, lags2, max2, peakInd2 = self.autoCorrCounts(self.dataA['dateTime'], 
                                    self.dataA['dos2rate'], 
                                    autoCorrRange, ax=self.ax[-2], 
                                    label='dos2rate', c='b')
                #self.corrCounts(self.dataA['dateTime'],self.dataB['dateTime'], 
#                                self.dataA['dos1rate'], #self.dataB['dos1rate'], 
#                                crossCorrRange, ax=self.ax[-1])
                dos12corr, dos12lags = self.corrCounts(
                                    self.dataA['dateTime'], 
                                    self.dataA['dateTime'], 
                                    self.dataA['dos1rate'], 
                                    self.dataA['dos2rate'], 
                            crossCorrRange,      ax=self.ax[-1])                            
                self.ax[-2].legend(loc=1)
                self.ax[-1].legend(loc=1)
                
                # First check that dos1 and dos12 are correlated, then
                # then if max(dos2) > 0.5*max(dos1) then it is noise. 
                if ( (len(np.where(peakInd1 == 4)[0]) == 1 or 
                        len(np.where(peakInd1 == 2)[0]) == 1 or 
                        len(np.where(peakInd2 == 4)[0]) == 1 or 
                        len(np.where(peakInd2 == 2)[0]) == 1) and
                        (max(self.dataA['dos2rate'][validIdA]) > 1000 or
                        max(dos12corr) >= 0.9) ):
                    self.noiseFlag = 'false'
                #elif max(dos12corr[len(dos12corr)//2-2 : len(dos12corr)//2+2]) > 0.9:
                #    self.noiseFlag = 'false'
                else:
                    self.noiseFlag = 'true'

                # Annotate the plot if the detection is a false positive
                if self.noiseFlag == 'false': 
                    self.ax[-1].text(0.5, 0.9, 'Flagged as Noise', va='top',
                        ha='right', transform=self.ax[-1].transAxes)
                
                # Print a True/False flag if the first peak is less than
                # a threshold.
                #lag4 = np.where(peakInd == 4)[0] # Check if there is a peak at 0.4 s lag.W
                #if len(lag4) == 1 and np.percentile(counts[validIdt], 95) > 1000:
                #    
                # Backwards to that 'false' means bad detection.
                #    self.noiseFlag = 'false' 
                #else:
                #    self.noiseFlag = 'true'

                saveDate = pltRange[0].replace(
                                microsecond=0).isoformat().replace(':', '')
                saveName = '{}_microburst_validation.png'.format(saveDate)
                
                plt.tight_layout()
                if self.uniqueFolders:
                    plt.savefig(os.path.join(self.saveDir, self.noiseFlag, saveName))
                else:
                    plt.savefig(os.path.join(self.saveDir, saveName))
                for a in self.ax:
                    a.clear()
        return
        
    def plotCounts(self, ax, timeA, timeB, doskey, tRange, pDet, sDet):
        """
        This method plots the narrow microburst time range from
        both spacecraft, and annotate with L, MLT, Lat, Lon.

        pDet and sDet are the primary and secondary spacececraft
        detection times.
        """
        sc_opt = ['A', 'B'] # Spacecraft options (for the plot labels.)
        dosA = self.dataA[doskey]
        dosB = self.dataB[doskey]

        validIdA = np.where((dosA != -1E31) & 
                            (timeA > tRange[0]) & 
                            (timeA < tRange[1]))[0]
        validIdB = np.where((dosB != -1E31) & 
                            (timeB > tRange[0]) & 
                            (timeB < tRange[1]))[0]
        if len(validIdA) == 0 or len(validIdB) == 0:
            return -1

        # Plot the timeseries.
        ax.plot(timeA[validIdA], dosA[validIdA], 'r', 
                    label='AC-6 {}'.format(self.primary_sc))
        ax.plot(timeB[validIdB], dosB[validIdB], 'b',
                    label='AC-6 {}'.format(sc_opt[sc_opt != self.primary_sc][0]))
        # Mark detection times
        #for m in pDet:
        #    ax.text(m, 0.1, 'x', color='r', transform=ax.transAxes)
        #for m in sDet:
        #    ax.text(m, 0.05, 'x', color='b', transform=ax.transAxes)
        maxCC = np.max([dosA[validIdA].max(), dosB[validIdB].max()])
        ax.scatter(pDet, maxCC*np.ones_like(pDet), c='r', marker='*', s=50) 
        ax.scatter(sDet, 0.9*maxCC*np.ones_like(sDet), c='b', marker='x', s=30)
        # for det in sDet:
        #     self.ax.axvline(det, c='b', lw=)
        ax.set(ylabel='{} [counts/s]'.format(doskey), xlabel='UTC',
            xlim=tRange,
            title='AC6-{} microburst validation | {}'.format(
                                self.primary_sc, tRange[0].date()))
        ax.locator_params(axis='x', nbins=3)

        #self.bx.set(xlabel='UTC', ylabel='Alpha [Deg] (dashed)')
        ax.legend(loc=1)

        meanFlag = np.mean(self.dataA['flag'][validIdA])
        textStr = ('L={} MLT={}\nlat={} lon={}\ndist={} LCT={}\nflag={}'.format(
            round(np.mean(self.dataA['Lm_OPQ'][validIdA])),
            round(np.mean(self.dataA['MLT_OPQ'][validIdA])),
            round(np.mean(self.dataA['lat'][validIdA])),
            round(np.mean(self.dataA['lon'][validIdA])),
            round(np.mean(self.dataA['Dist_In_Track'][validIdA])),
            round(np.mean(self.dataA['Loss_Cone_Type'][validIdA])),
            round(meanFlag)))
        ax.text(0.05, 0.95, textStr, transform=ax.transAxes, 
            va='top')
        return 1
        
    def autoCorrCounts(self, times, counts, tRange, ax=None, norm=True, label='dos1rate', c='r'):
        """
        This function calculates the autocorrelation of the counts array in
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
        peakInd, properties = scipy.signal.find_peaks(ac, prominence=0.1)
        
        if ax is not None:
            ax.plot(lags, ac, c=c, label='{} auto-correlation'.format(label))
            ax.set(xlabel='Lag [s]', ylabel='correlation coefficient', ylim=(-1, 1))
            ax.scatter(lags[peakInd], ac[peakInd], marker='+')

            # Mark prominances
            ax.vlines(x=lags[peakInd], ymax=ac[peakInd], 
                    ymin=ac[peakInd] - properties["prominences"], color="C1")
            for i, pi in enumerate(peakInd):
                ax.text(lags[pi], ac[pi] - properties["prominences"][i] - 0.2, 
                        round(properties["prominences"][i], 2), 
                        va='bottom', ha='center', color='C1')
        return ac, lags, max(counts[validIdt]), peakInd

    def corrCounts(self, timesA, timesB, countsA, countsB, tRange, ax=None, norm=True):
        """
        This method calculates the autocorrelation of the counts array in
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
        #ac = ac[ac.size//2:]
        # Lags are normalized to seconds  
        lags = np.arange(-cc.size/2, cc.size/2)/10 
        
        if norm:
            cc /= np.sqrt(len(x)*np.var(x)*len(y)*np.var(y))
            
        # Identify peaks
        peakInd, _ = scipy.signal.find_peaks(cc)
        
        if ax is not None:
            ax.plot(lags, cc, c='r', label='cross correlation')
            ax.axhline(0.7, c='g', ls=':', label='0.7')
            ax.axhline(0.8, c='g', ls='--', label='0.8')
            ax.axhline(0.9, c='g', ls='-', label='0.9')
            ax.set(xlabel='Lag [s] (shift blue trace by...)', ylabel='correlation coefficient', ylim=(-1, 1))
            # Mark peaks
            ax.scatter(lags[peakInd], cc[peakInd], marker='+') 
        return cc, lags
        
    def _load_catalog_(self, fPath):
        """
        This method reads in a catalog csv file, and saves it to a dictionary.
        """
        with open(fPath) as f:
            r = csv.reader(f)
            keys = next(r)
            rawData = np.array(list(r))
        data = {key:rawData[:, i] for i, key in enumerate(keys)}

        # Now convert the times array(s) to datetime, and all others to float
        for key in itertools.filterfalse(lambda x:'dateTime' in x, 
                                        data.keys()):
                data[key] = data[key].astype(float)
        for key in filter(lambda x: 'dateTime' in x, data.keys()):
            data[key] = np.array([dateutil.parser.parse(i) for 
                                i in data[key]])
        return data 
        
    def _findDateBounds(self):
        """ 
        This function will calculate the unique dates in the catalogue file
        to speed up the plotting reutine.
        """
        dates = np.array([t.date() for t in self.cDataPrimary['dateTime']])
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
    sc_id = 'A'
    catPathA = ('/home/mike/research/ac6-microburst-scale-sizes/'
                'data/microburst_catalogues/AC6A_microbursts_v3.txt')
    catPathB = ('/home/mike/research/ac6-microburst-scale-sizes/'
                'data/microburst_catalogues/AC6B_microbursts_v3.txt')
    p = ValidateDetections(sc_id, catPathA, catPathB)
    p.sep_range = [400, 500]
    p.plotLoop()
    p.sep_range = [0, 20]
    p.plotLoop()
