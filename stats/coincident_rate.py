import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
from datetime import datetime, timedelta
import dateutil.parser
import csv
import itertools
import scipy.signal
#import sys
import os 
import functools

#sys.path.insert(0, '/home/mike/research/mission_tools/ac6')
import mission_tools.ac6.read_ac_data as read_ac_data

class OccuranceRate:
    def __init__(self, sc_id, date, catV, catPath=None):
        """
        NAME: OccuranceRate
        USE:  Calculates the occurance rate of microbursts
              from a single spacecraft as a function of
              radiation belt pass. Radiation belt pass is
              defined as 4 < L < 8 (kwargs for the 
              radBeltIntervals method.) Radiation belt 
              interval indicies are stored in 
              self.intervals array(nPasses x 2). Occurance 
              rate for each pass is saved in self.rates
              array. Occurance rate is calculated by adding 
              up all of the microbursts widths for each pass
              (assuming a width function for which two are
              defined and are chosen by the mode kwarg of 
              the occurance_rate method) and then divided 
              by the pass time.
        INPUT: sc_id: spacecraft id
               date:  data to load data from
               catV:  microburst catalog version
               catPath: microburst catalog path (if not
                        specied, the __init__ method 
                        will use a default path)
        AUTHOR: Mykhaylo Shumko
        RETURNS: self.intervals: array(nPasses x 2) of 
                    radiation belt pass start/end indicies
                 self.rates: array(nPasses) of the 
                    microburst duty cycle between 0 and 1 
                    for each pass.
        MOD:     2018-10-23
        """
        self.sc_id = sc_id
        self.date = date
        self._load_sc_data()

        if catPath is None:
            catPath = ('/home/mike/research/ac6_microburst_scale_sizes'
                      '/data/microburst_catalogues')
        self._load_catalog(catPath, catV)
        return

    def radBeltIntervals(self, lbound=4, ubound=8, dawnFilter=False):
        """
        This method finds all of the radiation belt passes in a day 
        and separates it into start/end indicies.
        """
        # Find the start/end indicies of each rad belt pass.
        L = np.abs(self.data['Lm_OPQ'])

        # Filter by L and MLT if dawnFilter is True.
        if dawnFilter:
            MLT = self.data['MLT_OPQ']
            idL = np.where((L > lbound) & (L < ubound) & 
                            (MLT > 0) & (MLT < 12))[0]
        else:
            idL = np.where((L > lbound) & (L < ubound))[0]
        # Find breaks in the rad belt indicies to separate 
        # the passes.
        conv = np.convolve([1, -1], idL, mode = 'valid') - 1
        consecutiveFlag = np.where(conv != 0)[0] + 1
        startInd = np.insert(consecutiveFlag, 0, 0)
        endInd = np.insert(consecutiveFlag, 
                len(consecutiveFlag), len(idL)-1)
        # Save pass indicies to self.intervals.
        self.intervals = np.zeros((len(startInd), 2), dtype=int)
        for i, (i_s, i_e) in enumerate(zip(startInd, endInd)):
            self.intervals[i, :] = [idL[i_s], idL[i_e-1]]
        return 

    def occurance_rate(self, mode='static', **kwargs):
        """ 
        This method calculates the microburst occurance rates
        for each pass assuming a microburst width assumption.
        The width is either a fixed value in seconds or 
        another class method.
        """
        verbose = kwargs.get('verbose', False)
        testPlot = kwargs.get('testPlot', False)
        rel_height = kwargs.get('rel_height', 1)
        
        # Assume a static microburst width.
        if mode == 'static':
            w = kwargs.get('static_width', 0.5)
            width = functools.partial(self.static_width, w)
        # Assume a width that is calculated by 
        # scipy.signal.peak_width() after detrending.
        elif mode == 'prominence':
            width = functools.partial(self.prominence_width, 
                                      rel_height=rel_height, 
                                      testPlot=testPlot)
            
        # Convert the catalog datetimes to numeric to enable
        # faster time comparison.
        nDetTime = date2num(self.cat['dateTime'])
        self.rates = np.zeros(self.intervals.shape[0], 
                                dtype=float)
        if testPlot:
            _, self.ax = plt.subplots(2, sharex=True)
        # Loop over the passes.
        for i, (i_start, i_end) in enumerate(self.intervals):                
            # Find microbursts in ith pass
            startTime = self.data['dateTime'][i_start]
            endTime = self.data['dateTime'][i_end]
            idt = np.where((nDetTime > date2num(startTime)) & 
                            (nDetTime < date2num(endTime)))[0]       
            if verbose: 
                print('Processing pass', i, 'len(idt)=', len(idt))

            # determine the duration of each microburst in the
            # current pass
            for t_i in idt:
                self.rates[i] += width(self.cat['dateTime'][t_i])
            pass_duration = (endTime - startTime).total_seconds()
            if verbose:
                print('rates =', self.rates[i], 
                    'duration=', pass_duration)
            self.rates[i] /= pass_duration
        return
        
    def static_width(self, width, t):
        """ Returns a constant width """
        return width
        
    def prominence_width(self, t, testPlot=True, rel_height=0.5,
                         window_width=1):
        """ 
        This method implements an algorithm to calculate the microburst 
        width based on its topological prominence.

        scipy.signal.find_peaks()
        """
        # Clip the data to the microburst time (with some window 
        # around it).
        lt = t - timedelta(seconds=window_width)
        ut = t + timedelta(seconds=window_width)
        indicies = np.where((self.data['dateTime'] > lt) & 
                           (self.data['dateTime'] < ut))[0]
        
        # detrend microburst time series
        iPeak = np.where(self.data['dateTime'][indicies] == t)[0]
        if len(iPeak) != 1: raise ValueError('0 or > 1 peaks found!')
        detrended = scipy.signal.detrend(self.data['dos1rate'][indicies])

        # Try to find the width of the center peak.
        try:
            width = scipy.signal.peak_widths(detrended, iPeak, 
                                            rel_height=0.5)
        # If it fails, search for nearest peak to the center time 
        # and try again.
        except ValueError as err:
            if 'is not a valid peak' in str(err):
                peaks, _ = scipy.signal.find_peaks(detrended, 
                                            prominence=None)
                iPeak = peaks[np.argmin(np.abs(peaks - len(indicies)//2))]
                width = scipy.signal.peak_widths(detrended, [iPeak],
                                                 rel_height=rel_height)
        # Make a test plot for debugging purposes.
        if testPlot:
            saveDir = '/home/mike/temp_plots'
            if not os.path.exists(saveDir):
                print('made dir', saveDir)
                os.makedirs(saveDir)
            self.ax[0].plot(self.data['dos1rate'][indicies])
            self.ax[0].plot(self.data['dos1rate'][indicies]-detrended)
            self.ax[1].plot(detrended)
            self.ax[1].hlines(*width[1:])
            self.ax[1].set_ylim(top=1.5*detrended[iPeak])
                
            plt.savefig(os.path.join(saveDir,
                        '{}_microburst_validation.png'.format(
                        t.replace(microsecond=0).isoformat()))
                        )
            for a in self.ax:
                a.cla()
        return 0.1*width[0][0]

    def _load_sc_data(self):
        """ Loads AC-6 10 Hz data """
        print('Loading AC6-{} data from {}'.format(
            self.sc_id.upper(), self.date.date()))
        self.data = read_ac_data.read_ac_data_wrapper(
                                self.sc_id, self.date)
        return

    def _load_catalog(self, catPath, v):
        """ 
        Loads the microburst catalog of version v spacecraft given 
        by self.sc_id. 
        """
        fPath = os.path.join(catPath, 
            'AC6{}_microbursts_v{}.txt'.format(self.sc_id.upper(), v))
        print('Loading catalog,', 'AC6{}_microbursts_v{}.txt'.format(
            self.sc_id.upper(), v))
        with open(fPath) as f:
            r = csv.reader(f)
            keys = next(r)
            rawData = np.array(list(r))
        self.cat = {}
        for (i, key) in enumerate(keys):
            self.cat[key] = rawData[:, i]

        # Now convert the times array(s) to datetime, 
        # and all others except 'burstType' to a float.
        timeKeys = itertools.filterfalse(
            lambda x:'dateTime' in x or 'burstType' in x, self.cat.keys())
        for key in timeKeys:
                self.cat[key] = self.cat[key].astype(float)
        for key in filter(lambda x: 'dateTime' in x, self.cat.keys()):
            self.cat[key] = np.array([dateutil.parser.parse(i) 
                for i in self.cat[key]])
        return

class CoincidenceRate:
    def __init__(self, date, catV, catPath=None):
        """
        NAME: CoincidenceRate
        USE:  Calculates the coincidence rate of microbursts
              from both spacecraft as a function of
              radiation belt pass. This classe uses the 
              methods from OccuranceRate class, so refer to
              those docs on how the occurance rate and 
              radiation belt passes are calculated.
        INPUT: date:  data to load data from
               catV:  microburst catalog version
               catPath: microburst catalog path (if not
                        specied, the __init__ method 
                        will use a default path)
        AUTHOR: Mykhaylo Shumko
        RETURNS: self.intervals: array(nPasses x X) of 
                    radiation belt pass start/end indicies
                    for both spacecraft
                 self.rates: array(nPasses) of the coincident 
                    microburst duty cycle between 0 and 1 
                    for each pass.
        MOD:     2018-10-23
        """
        # Load AC6-A data and microburst catalog
        self.occurA = OccuranceRate('A', date, catV, catPath=catPath)
        # Load AC6-B data and microburst catalog.
        self.occurB = OccuranceRate('B', date, catV, catPath=catPath)
        return
    
    def radBeltIntervals(self):
        """
        NAME:  radBeltIntervals
        USE:   This method calculates the radiation belt 
               passes from both spacecraft.
        INPUT: None
        AUTHOR: Mykhaylo Shumko
        RETURNS: self.passes - an array(nPasses x 4) where
                            the columns indicate, in order,
                            start/stop indicies for AC6A, 
                            followed by indicies for AC6B 
                            for the matching passes.
        MOD:     2018-10-25
        """
        # Calculate rad belt pass intervals from each 
        # spacecraft.
        self.occurA.radBeltIntervals()
        self.occurB.radBeltIntervals()

        self.passes = np.zeros((0, 4), dtype=object)

        tAn = date2num(cr.occurA.data['dateTime'][cr.occurA.intervals])
        tBn = date2num(cr.occurB.data['dateTime'][cr.occurB.intervals])

        # Find matching indicies by looping over start times of AC6A 
        # pass interval array
        for i, t in enumerate(tAn[:, 0]):
            # For each start time in AC6A, look for a corresponding 
            # start times in AC6B data within 5 minutes
            dt = np.abs(tBn[:, 0] - t)
            idmin = np.argmin(dt)
            if dt[idmin] < sec2day(5*60):
                newRow = np.concatenate((self.occurA.intervals[i], 
                                        self.occurB.intervals[idmin]
                                        ))
                self.passes = np.vstack((self.passes, newRow))
        return

    def sortBursts(self, ccThresh=0.8, ccOverlap=2):
        """
        NAME:   sortBursts
        USE:    This method loops through every pass for which there
                is data from both spacecraft, and for each detecton
                made by either AC6A or AC6B, checks if there is a 
                correlated peak at the same time, and in the same
                position. Uses cross-correlations (CC)
        INPUT:  ccThresh = 0.8 - CC threshold for detection.
                ccOverLap = 1  - CC overlap, to account for Poisson 
                                 statistics for short-duration 
                                 microbursts and small timing offsets
                                (little to none expected).
        AUTHOR: Mykhaylo Shumko
        RETURNS: self.tBurst - array(nDet, 2) of temporal microbursts
                 self.sBurst - array(nDet, 3) of curtains with AC6A 
                               time, AC6B time, and CC.
                 self.aBurst - array(nDet, 2) of ambigious bursts
                 First column is time, second is cross correlation
                 value.
        MOD:     2018-10-25
        """
        self.tBurst = np.zeros((0, 2), dtype=object) 
        self.sBurst = np.zeros((0, 4), dtype=object) 
        self.aBurst = np.zeros((0, 2), dtype=object) 
        
        # Convert detection times to numbers for easy comparison.
        catAtimes = date2num(self.occurA.cat['dateTime'])
        catBtimes = date2num(self.occurB.cat['dateTime'])

        # Loop over passes
        for i, (passAi, passAf, passBi, passBf) in enumerate(self.passes):
            print('Processing pass', i)
            # Find ime range for this pass. 
            aRange = date2num([self.occurA.data['dateTime'][passAi],
                             self.occurA.data['dateTime'][passAf]])
            bRange = date2num([self.occurB.data['dateTime'][passBi],
                             self.occurB.data['dateTime'][passBf]])
            # Find all bursts separately detected on AC6A/B for 
            # this pass
            burstsA = np.where((catAtimes > aRange[0]) & 
                            (catAtimes < aRange[1]))[0]
            burstsB = np.where((catBtimes > bRange[0]) & 
                            (catBtimes < bRange[1]))[0]
            print('Found {} AC6A bursts and {} AC6B bursts'.format(
                len(burstsA), len(burstsB)))

            # Loop over bursts from AC6A.              
            for bA in burstsA:
                flag, cc, t2, timeLag  = self.cross_correlate('A', bA, ccOverlap)
                # Save to the appropriate array. 
                if flag == 't':
                    self.tBurst = np.vstack((self.tBurst, 
                    [self.occurA.cat['dateTime'][bA], cc]))
                # Likeliy need to do something special 
                # with curtains. 
                elif flag == 's':
                    self.sBurst = np.vstack((self.sBurst, 
                    [self.occurA.cat['dateTime'][bA], t2, timeLag, cc]))
                elif flag == 'a':
                    self.aBurst = np.vstack((self.aBurst, 
                    [self.occurA.cat['dateTime'][bA], cc]))

            # Loop over bursts from AC6B.              
            for bB in burstsB:
                flag, cc, t2, timeLag = self.cross_correlate('B', bB, ccOverlap)
                # Save to the appropriate array. 
                if flag == 't':
                    self.tBurst = np.vstack((self.tBurst, 
                    [self.occurB.cat['dateTime'][bB], cc]))
                elif flag == 's':
                    self.sBurst = np.vstack((self.sBurst, 
                    [t2, self.occurB.cat['dateTime'][bB], timeLag, cc]))
                elif flag == 'a':
                    self.aBurst = np.vstack((self.aBurst, 
                    [self.occurB.cat['dateTime'][bB], cc]))
        
        return
        
    def sortArrays(self):
        """ 
        This method sorts and removes duplicated from the sorted
        burst arrays and.
        """
        # Remove duplicates. This only removes elements whenever
        # both time and cc are the same (true in this case since
        # CC is symmetric).
        self.tBursts = list({ (i[0], i[1]) for i 
                        in self.tBursts})
        self.sBursts = list({ (i[0], i[1]) for i 
                        in self.sBursts})
        self.aBursts = list({ (i[0], i[1]) for i 
                        in self.aBursts})
        
        # Sort times
        self.tBurst = np.array(sorted(
                               self.tBurst, key=lambda x : x[0]))
        self.sBurst = np.array(sorted(
                               self.sBurst, key=lambda x : x[0]))
        self.aBurst = np.array(sorted(
                               self.aBurst, key=lambda x : x[0]))
        return 

    def cross_correlate(self, sc_id, i, ccOverlap, ccWindow=0.5):
        """
        NAME:   cross_correlate
        USE:    This method takes in a microburst indicie, and
                calculates the cross-correlation for that time 
                series with the time series of the other 
                spacecraft at the same time, and same position.
                
        INPUT:  sc_id     - Spacecraft id. Either 'A' or 'B'.
                i         - Index of the detection for sc_id.
                ccOverlap - Cross correlation overlap.
                ccWindow  - Cross correlation window width.
        AUTHOR: Mykhaylo Shumko
        RETURNS: flag - Detection type. Can be 't' for temporal
                        peak, 's' for spatial peak, and 'a' for 
                        ambigrious.
                 cc   - Max cross-correlation value
        MOD:     2018-10-25
        """
        # Time cross-correlation coefficient
        ccT = self._correlate_time(sc_id, i, ccOverlap, ccWindow)
        # Space cross-correlation coefficient
        t, t2, ccS, timeLag = self._correlate_space(sc_id, i, ccOverlap, ccWindow)
        # Logic on if the detection is a microburst, curtain, or 
        # ambigious 
        if ccT == ccS:
            return 'a', ccT, None, None
        elif ccT > ccS:
            return 't', ccT, None, None
        else:
            return 's', ccS, t2, timeLag

    def _correlate_time(self, sc_id, i, ccOverlap, ccWindow=2):
        """ Cross correlate data at the same times """
        # Find the center time of that burst
        if sc_id.upper() == 'A':
            t = self.occurA.cat['dateTime'][i]
        else:
            t = self.occurB.cat['dateTime'][i]
        # Cross correlate in time.
        tWindowA = np.where(
                (cr.occurA.data['dateTime'] > t-timedelta(seconds=ccWindow/2)) &  
                (cr.occurA.data['dateTime'] < t+timedelta(seconds=ccWindow/2)) &
                (cr.occurA.data['dos1rate'] != -1E31)
                )[0]     
        tWindowB = np.where(
                (cr.occurB.data['dateTime'] > t-timedelta(seconds=ccWindow/2)) &  
                (cr.occurB.data['dateTime'] < t+timedelta(seconds=ccWindow/2)) &
                (cr.occurB.data['dos1rate'] != -1E31)
                )[0]
        # Widen the CC window.
        tWindowB = np.insert(tWindowB, 0, tWindowB[0]-ccOverlap//2)
        tWindowB = np.append(tWindowB, tWindowB[0]+ccOverlap//2)

        x = cr.occurA.data['dos1rate'][tWindowA] - cr.occurA.data['dos1rate'][tWindowA].mean()
        y = cr.occurB.data['dos1rate'][tWindowB] - cr.occurB.data['dos1rate'][tWindowB].mean()
        ccArr = np.correlate(x, y, mode='same')/np.sqrt(len(x)*len(y)*np.var(x)*np.var(y))
        return max(ccArr)

    def _correlate_space(self, sc_id, i, ccOverlap, ccWindow=0.5):
        """
        Cross correlate data at the same place.
        """
        # Depending on the spacecraft, find the correct indicies 
        # for the spatial CC.
        if sc_id.upper() == 'A':
            t = self.occurA.cat['dateTime'][i] # A center time
            timeLag = self.occurA.cat['Lag_In_Track'][i]
            t2 = t - timedelta(seconds=timeLag) # B center time.
            tWindowA = np.where(
                (cr.occurA.data['dateTime'] > t-timedelta(seconds=ccWindow/2)) &  
                (cr.occurA.data['dateTime'] < t+timedelta(seconds=ccWindow/2)) &
                (cr.occurA.data['dos1rate'] != -1E31)
                )[0]     
            tWindowB = np.where(
                (cr.occurB.data['dateTime'] > t-timedelta(seconds=ccWindow/2-timeLag)) &  
                (cr.occurB.data['dateTime'] < t+timedelta(seconds=ccWindow/2-timeLag)) &
                (cr.occurB.data['dos1rate'] != -1E31)
                )[0]
        else:
            t = self.occurB.cat['dateTime'][i] # B center time
            timeLag = self.occurB.cat['Lag_In_Track'][i]
            t2 = t + timedelta(seconds=timeLag) # A center time.
            tWindowA = np.where(
                (cr.occurA.data['dateTime'] > t-timedelta(seconds=ccWindow/2+timeLag)) &  
                (cr.occurA.data['dateTime'] < t+timedelta(seconds=ccWindow/2+timeLag)) &
                (cr.occurA.data['dos1rate'] != -1E31)
                )[0]     
            tWindowB = np.where(
                    (cr.occurB.data['dateTime'] > t-timedelta(seconds=ccWindow/2)) &  
                    (cr.occurB.data['dateTime'] < t+timedelta(seconds=ccWindow/2)) &
                    (cr.occurB.data['dos1rate'] != -1E31)
                    )[0]

        # Widen the CC window.
        tWindowB = np.insert(tWindowB, 0, tWindowB[0]-ccOverlap//2)
        tWindowB = np.append(tWindowB, tWindowB[0]+ccOverlap//2)

        # Correlate the two signals
        x = cr.occurA.data['dos1rate'][tWindowA] - cr.occurA.data['dos1rate'][tWindowA].mean()
        y = cr.occurB.data['dos1rate'][tWindowB] - cr.occurB.data['dos1rate'][tWindowB].mean()
        ccArr = np.correlate(x, y, mode='same')/np.sqrt(len(x)*len(y)*np.var(x)*np.var(y))
        return t, t2, max(ccArr), timeLag

    def test_plots(self, window=5, saveDir='/home/mike/temp_plots'):
        """ 
        This method plots detections and their temporal and spatial
        cross-correlations.
        """
        from matplotlib.backends.backend_pdf import PdfPages
        _, ax = plt.subplots(2, sharex=True)  

        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
            print('Made plot directory:', saveDir)

        # Iterate over the microburst detections.
        savePath = os.path.join(saveDir, '{}_microburst_validation.pdf'.format(
                    datetime.now().date()))
        with PdfPages(savePath) as pdf:
            for (t, cc) in self.tBurst:
                self._test_plots_microbursts(ax, t, cc, window)
                #self._test_plots_space(ax[1], tA, tB, lag, cc, window)
                pdf.savefig()    
                for a in ax:
                    a.cla()
        # # Iterate over the curtain detections.            
        # savePath = os.path.join(saveDir, '{}_curtain_validation.pdf'.format(
        #             datetime.now().date()))
        # with PdfPages(savePath) as pdf:
        #     for (tA, tB, lag, cc) in self.sBurst:
        #         self._test_plots_space(ax[1], tA, tB, lag, cc, window)
        #         pdf.savefig()    
        #         for a in ax:
        #             a.cla()
        return

    def _test_plots_microbursts(self, ax, t, cc, window):
        """ 
        This helper method plots the microburst time series with 
        the shifted and unshifted data.
        """
        idtA = np.where(
            (self.occurA.data['dateTime'] > t-timedelta(seconds=window/2)) &  
            (self.occurA.data['dateTime'] < t+timedelta(seconds=window/2)) &
            (self.occurA.data['dos1rate'] != -1E31)
            )[0]     
        idtB = np.where(
            (self.occurB.data['dateTime'] > t-timedelta(seconds=window/2)) &  
            (self.occurB.data['dateTime'] < t+timedelta(seconds=window/2)) &
            (self.occurB.data['dos1rate'] != -1E31)
            )[0]
        dt = self.occurB.data['Lag_In_Track'][idtB[0]] # in-track lag
        # Shift AC6B times by dt.
        bTimes = np.array([t - dt for t in self.occurB.data['dateTime'][idtB]])

        idtB_shifted = np.where(
            (self.occurB.data['dateTime'] > t-timedelta(seconds=window/2)) &  
            (self.occurB.data['dateTime'] < t+timedelta(seconds=window/2)) &
            (self.occurB.data['dos1rate'] != -1E31)
            )[0]

        ax[0].plot(self.occurA.data['dateTime'][idtA], 
                    self.occurA.data['dos1rate'][idtA], 
                    'r', label='AC6A')
        ax[0].plot(self.occurB.data['dateTime'][idtB], 
                    self.occurB.data['dos1rate'][idtB], 
                    'b', label='AC6B')
        ax[0].text(0.1, 0.9, 'CC={}'.format(cc), 
                    transform=ax[0].transAxes)
        ax[0].set(title='{} | AC6 microburst validation'
                    ''.format(t.replace(microsecond=0).isoformat()))
        ax[0].axvline(t)
        ax[0].legend(loc=1)
        return

    # def _test_plots_space(self, ax, tA, tB, lag, cc, window):
    #     """ This helper method plots the time-aligned time series."""
    #     idtA = np.where(
    #         (self.occurA.data['dateTime'] > t-timedelta(seconds=window/2)) &  
    #         (self.occurA.data['dateTime'] < t+timedelta(seconds=window/2)) &
    #         (self.occurA.data['dos1rate'] != -1E31)
    #         )[0]     
    #     idtB = np.where(
    #         (self.occurB.data['dateTime'] > t-timedelta(seconds=window/2)) &  
    #         (self.occurB.data['dateTime'] < t+timedelta(seconds=window/2)) &
    #         (self.occurB.data['dos1rate'] != -1E31)
    #         )[0]

    #     ax.plot(self.occurA.data['dateTime'][idtA], 
    #                 self.occurA.data['dos1rate'][idtA], 
    #                 label='AC6A')
    #     ax.plot(self.occurB.data['dateTime'][idtB], 
    #                 self.occurB.data['dos1rate'][idtB], 
    #                 label='AC6B')
    #     ax.text(0.1, 0.9, 'CC={}'.format(cc), 
    #                 transform=ax[0].transAxes)
    #     ax.set(title='{} | AC6 microburst validation'
    #                 ''.format(t.replace(microsecond=0).isoformat()))
    #     ax.axvline(t)
    #     ax.legend(loc=1)
    #     return

def sec2day(s):
    """ Convert seconds to fraction of a day."""
    return s/86400

if __name__ == '__main__':
    cr = CoincidenceRate(datetime(2016, 10, 14), 3)
    #cr = CoincidenceRate(datetime(2015, 8, 28), 3)
    cr.radBeltIntervals()
    cr.sortBursts()
    cr.test_plots()

    # plt.scatter(cr.tBurst[:, 0], cr.tBurst[:, 1])
    # plt.xlim(cr.tBurst[0, 0], cr.tBurst[-1, 0]); plt.show()
    #print(cr.occrA.intervals)
    #print(cr.occrB.intervals)

    # tA = np.array([t.isoformat() for t in 
    #         cr.occurA.data['dateTime'][cr.occurA.intervals]])
    # tB = np.array([t.isoformat() for t in 
    #         cr.occurB.data['dateTime'][cr.occurB.intervals]])
    # np.savetxt('ac6a_pass_times.csv', tA, delimiter=',')
    # np.savetxt('ac6b_pass_times.csv', tB, delimiter=',')

    # with open('ac6a_pass_times.csv', 'w') as f:
    #     w = csv.writer(f)
    #     w.writerows(cr.occurA.data['dateTime'][cr.occurA.intervals])

    # with open('ac6b_pass_times.csv', 'w') as f:
    #     w = csv.writer(f)
    #     w.writerows(cr.occurB.data['dateTime'][cr.occurB.intervals])

    # ### CODE TO TEST THE OCCURANCE RATE TIMESERIES ###
    # o = OccuranceRate('A', datetime(2016, 10, 14), 3)
    # o.radBeltIntervals()
    # #for h in np.arange(0.5, 1.1, 0.5):
    # o.occurance_rate(mode='prominence')
    # #print('rel_height =', h)
    
    # ### PLOT OCCURANCE RATE TIMESERIES ###
    # startTimes = o.data['dateTime'][o.intervals[:, 0]]
    # _, ax = plt.subplots(2, sharex=True)
    # ax[0].scatter(startTimes, o.rates)
    # ax[0].set_ylabel('Microbust duty cycle\n(microburst/pass time)')
    
    # idt = np.where((o.cat['dateTime'] > startTimes[0]) & 
    #                 (o.cat['dateTime'] < startTimes[-1]))[0]

    # ax[1].scatter(o.cat['dateTime'][idt], o.cat['AE'][idt])

    # ax[1].set_ylabel('AE (nT)')
    # ax[-1].set_xlabel('UTC')
    # plt.tight_layout()
    # plt.show()

    # tAn = date2num(cr.occurA.data['dateTime'][cr.occurA.intervals])
    # tBn = date2num(cr.occurB.data['dateTime'][cr.occurB.intervals])
    # tA = cr.occurA.data['dateTime'][cr.occurA.intervals]
    # tB = cr.occurB.data['dateTime'][cr.occurB.intervals]

    # # Find matching indicies
    # # Loop pver start times of one array
    # for i, t in enumerate(tAn[:, 0]):
    #     # Look for start times of second array within 5 minutes
    #     dt = np.abs(tBn[:, 0] - t)
    #     idmin = np.argmin(dt)
    #     if dt[idmin] < sec2day(5*60):
    #         print('Found match! (i, t_A) =', i, tA[i, 0], 
    #             '(idmin, tB[idmin]) =', idmin, tB[idmin, 0] )
    #     else:
    #         print('No match (i, t_A) =', i, tA[i, 0], 
    #             '(idmin, tB[idmin]) =', idmin, tB[idmin, 0] )

