import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
from datetime import datetime, timedelta
import dateutil.parser
import csv
import itertools
import scipy.signal
import pandas as pd
import os 
import functools

#sys.path.insert(0, '/home/mike/research/mission_tools/ac6')
import mission_tools.ac6.read_ac_data as read_ac_data

class CumulativeDist:
    def __init__(self, catV, catPath=None, cc_width=1, cc_overlap=0.2, verbose=True):
        """
        NAME: CumulativeDist
        USE:  For each day, the loop method calculates the 
              cross-correlation for spatial and temporal data.
              This class is inherited from CoincidenceRate and
              the main method that is overwritten is loop()
        INPUT: sc_id: spacecraft id
               date:  data to load data from
               catV:  microburst catalog version
               catPath: microburst catalog path (if not
                        specied, the __init__ method 
                        will use a default path)
        AUTHOR: Mykhaylo Shumko
        RETURNS: Appends the temporal and spatial cross-correlations 
                to the catalog file. 
        MOD:     2019-01-06
        """
        if catPath is None:
            catPath = ('/home/mike/research/ac6_microburst_scale_sizes'
                      '/data/microburst_catalogues')
        self.catA = self._load_catalog('A', catPath, catV)
        self.catB = self._load_catalog('B', catPath, catV)
        self.catV = catV
        self.cc_width = cc_width
        self.cc_overlap=cc_overlap
        self.verbose = verbose
        return

    def loop(self):
        """ 
        Loop over every detection and cross-correlate time and space. 

        LOOP OUTLINE:
        1. Loop over every date for which microbursts were observed at
            one time or another by both units.
        2. Run two separate loops to find all of the miocrobursts 
            observed in AC6A and AC6B catalogs for that day.
        3. For each microburst, find the indicies of 10 Hz data
            at the same time and position. If no data exists at 
            the same time or space, skip this event. Otherwise
            cross-correlate in space and time with a window width 
            given by cc_width.
        4. Save the cross-correlation values to self.data array.
        """
        # self.data contains a copy of the catalog values as well as
        # four values: time cross-correlation, space cross-correlation
        # center spatial time for unit A, center spatial time for unit B.
        # To identify which spacecraft the microburst was detected, look
        # for the spatial time which matches the 'dateTime' time.
        self.data = np.nan*np.zeros((0, len(c.catA.keys())+6))
        self.dates = self._find_loop_dates()

        for date in self.dates:
            if self.verbose: print('Processing data on', date.date())
            # Load 10 Hz data from both spacecraft.
            self.tenHzA = read_ac_data.read_ac_data_wrapper(
                                'A', date)
            self.tenHzB = read_ac_data.read_ac_data_wrapper(
                                'B', date)
            # Find the indicies of the microburst detections on date.
            idA = np.where(date.date() == self.catDatesA)[0]
            idB = np.where(date.date() == self.catDatesB)[0]
            # Loop over the daily detections and cross-correlate.
            if self.verbose: print('Looping over AC6-A detections')
            self._daily_microburst_loop('A', idA)
            if self.verbose: print('Looping over AC6-B detections')
            self._daily_microburst_loop('B', idB)            
        return

    def save_catalog(self, savePath=None):
        """ 
        Saves the catalog to a savePath file which has a default 
        name if not specified.
        """
        if savePath is None:
            saveDir = './../data/coincident_microbursts_catalogues'
            saveName = 'AC6_coincident_microbursts_v{}.txt'.format(self.catV)
            savePath = os.path.join(saveDir, saveName)
        columns = np.concatenate((list(self.catA.keys()), 
                    ['time_cc', 'space_cc', 
                    'time_spatial_A', 'time_spatial_B', 
                    'peak_width_A', 'peak_width_B']
                    ))
        df = pd.DataFrame(data=self.data, columns=columns)    
        # Add column with truncated seconds to help remove duplicates
        df['dateTime_trunc'] = df['dateTime'].apply(
                                lambda t: t.replace(microsecond=0))
        # Sort the truncated dateTime values (needed since for each day
        # AC6-A detections get added first and then AC6-B detections).
        df = df.sort_values(by='dateTime_trunc')
        # Now remove duplicate rows
        df = df.drop_duplicates('dateTime_trunc')
        df.drop('dateTime_trunc', axis=1, inplace=True)

        df.to_csv(savePath, index=False, na_rep='nan')
        return

    def _find_loop_dates(self):
        """ Get all of the unique AC6 data dates. """
        self.catDatesA = np.array([ti.date() 
                                for ti in self.catA['dateTime']])
        self.catDatesB = np.array([ti.date() 
                                for ti in self.catB['dateTime']])
        self.catNumDatesA = date2num(self.catDatesA)
        self.catNumDatesB = date2num(self.catDatesB)
        unique_num_dates = sorted(list(set(self.catNumDatesA) & 
                                        set(self.catNumDatesB)))
        return num2date(unique_num_dates)

    def _sort_catalog(self, duplicate_thresh=0.3):
        """ 
        This method sorts the catalog file to be saved and removes
        duplicate events that are defined as events occuring within a
        duplicate_thresh (in seconds)
        """
        data_sorted = sorted(self.data, key=lambda x:x[0])
        t = np.array([row[0] for row in data_sorted])
        dt = np.array([i.total_seconds() 
                    for i in (t[1:] - t[:-1])])
        idt = np.where(dt > duplicate_thresh)[0] # Non-duplicate values
        self.data = self.data[idt, :]
        return

    def _daily_microburst_loop(self, sc_id, idx):
        """
        Loops over the daily microbursts from a 
        spacecraft sc_id and over indicies idx. 
        For each microburst, the timeseries is 
        cross-correlated in space and time.
        """
        for i in idx:
            
            # If the Lag_In_Track value is an error, 
            # ignore that iteration.
            if sc_id.upper() == 'A':
                if np.isnan(self.catA['Lag_In_Track'][i]):
                    continue
            else:
                if np.isnan(self.catB['Lag_In_Track'][i]):
                    continue
            # Find time and space indicies.
            # out tuple consists of idtA, idtB, idtA_shifted, 
            # idtB_shifted, t0, t_sA, t_sB
            out = self._get_time_space_indicies(sc_id, i)
            
            # Check if time and space indicies are not empty.
            if len(out[0]) == 0 or len(out[1]) == 0 or len(out[2]) == 0 or len(out[3]) == 0:
                continue
            # If there is enough indicies, now cross-correlate
            time_cc = self.CC(out[0], out[1])
            space_cc = self.CC(out[2], out[3])
            peak_width_A, peak_width_B = self._find_peaks(out[0], out[1])
            # Need to think twice about the width in the spatial domain.
            #width_A, width_B = self._find_peaks(out[2], out[3])
            if sc_id.upper() == 'A':
                line = np.concatenate(([self.catA[key][i] for key in self.catA.keys()],
                                       [time_cc, space_cc, out[-2], out[-1], 
                                       peak_width_A, peak_width_B]))
            else:
                line = np.concatenate(([self.catB[key][i] for key in self.catB.keys()],
                                       [time_cc, space_cc, out[-2], out[-1], 
                                       peak_width_A, peak_width_B]))
            self.data = np.vstack((self.data, line))
        return

    def CC(self, iA, iB):
        """ 
        This method calculates the normalized Pearson's cross-correlation 
        coefficient between two AC6 time series indexed by iA and iB.
        """
        norm = np.sqrt(len(self.tenHzA['dos1rate'][iA])*\
                       len(self.tenHzB['dos1rate'][iB])*\
                       np.var(self.tenHzA['dos1rate'][iA])*\
                       np.var(self.tenHzB['dos1rate'][iB])) 
        # Mean subtraction.
        x = (self.tenHzA['dos1rate'][iA] - 
            self.tenHzA['dos1rate'][iA].mean() )
        y = (self.tenHzB['dos1rate'][iB] - 
            self.tenHzB['dos1rate'][iB].mean() )
        # Cross-correlate
        ccArr = np.correlate(x, y, mode='valid')
        # Normalization
        ccArr /= norm
        return max(ccArr)

    def _find_peaks(self, iA, iB, sample_thresh=1, peak_kwargs={}, smooth=3):
        """
        This method calls scipy.signal.find_peaks to attempt to
        find a peak within time_thresh of the center of the index
        array, iA and iB for indices from sc A and B, respectively.
        """
        if smooth > 1:
            countsA = np.convolve(np.ones(smooth)/smooth, self.tenHzA['dos1rate'][iA], 
                mode='same')
            countsB = np.convolve(np.ones(smooth)/smooth, self.tenHzB['dos1rate'][iB], 
                mode='same')
        # Find the peaks
        peaksA, propertiesA = scipy.signal.find_peaks(countsA, 
                **peak_kwargs, width=(None, None))
        peaksB, propertiesB = scipy.signal.find_peaks(countsB, 
                **peak_kwargs, width=(None, None))
        # The template sets just contain the indicies of iA and iB that 
        # are at the center of iA and iB +/- sample_thresh
        peak_template_A = set(np.arange(
                            len(iA)//2-sample_thresh, len(iA)//2+sample_thresh+1
                            ))
        peak_template_B = set(np.arange(
                            len(iB)//2-sample_thresh, len(iB)//2+sample_thresh+1
                            ))
        # Now find if there is a peak at or near the center.
        peak_df_A = pd.DataFrame({'ipeak':peaksA, 'width':propertiesA['widths']})
        peak_df_B = pd.DataFrame({'ipeak':peaksB, 'width':propertiesB['widths']})
        valid_peaks_A = peak_df_A[peak_df_A['ipeak'].isin(peak_template_A)]
        valid_peaks_B = peak_df_B[peak_df_B['ipeak'].isin(peak_template_B)]
        # If only one peak was found near the center, return the width_A/B. If 
        # none or more than 1 peaks were found, then return np.nan.
        if len(valid_peaks_A) == 1:
            width_A = valid_peaks_A.iloc[0, 1]/10
        else:
            width_A = np.nan
        if len(valid_peaks_B) == 1:
            width_B = valid_peaks_B.iloc[0, 1]/10
        else:
            width_B = np.nan
        return width_A, width_B

    def _get_time_space_indicies(self, sc_id, i):
        """ 
        This method calculates the time aligned and space-aligned
        indicies for cross-correlation. As for the time shifts,
        refer to AC6 data README which specifies that the 
        Lag_In_Track == t_A - t_B.

        OUTPUT:
            first two outputs - time-aligned indicies for 
                                AC6A and AC6B, respectfully.
            second two outputs - space-aligned indicies for 
                                AC6A and AC6B, respectfully.
            last three outputs - time-aligned time, shifted 
                                time for AC6A (same as 
                                time-aligned time if sc_id == 'A')
                                and last output is shifted time
                                for AC6B which will be the shifted 
                                center time if sc_id == A.
        """
        dt = timedelta(seconds=self.cc_width/2)
        overlapW = timedelta(seconds=self.cc_overlap/2)
        # Find the correct center time from either AC6A or B.
        if sc_id.upper() == 'A':
            t0 = self.catA['dateTime'][i]
        else:
            t0 = self.catB['dateTime'][i]

        #### Get temporally-aligned indicies ######
        idtA = np.where(  (self.tenHzA['dateTime'] >= t0-dt) 
                        & (self.tenHzA['dateTime'] < t0+dt)
                        & (self.tenHzA['dos1rate'] != -1E31))[0]
        # The overlapW windens the index array to accomidate a CC lag.
        idtB = np.where(  (self.tenHzB['dateTime'] >= t0-dt-overlapW) 
                        & (self.tenHzB['dateTime'] < t0+dt+overlapW)
                        & (self.tenHzB['dos1rate'] != -1E31))[0]

        ####### Get spatially aligned indicies ######
        # If sc_id is 'A', then, the spatial indicies are the 
        # same for AC6A, but shifted for AC6B. Since 
        #               lag = t_A - t_B
        # this implies that t_B_shifted = t_A - lag. If sc_id 
        # is 'B', then everything is opposite, i.e. 
        # t_A_shifted = t_B + lag
        if sc_id.upper() == 'A':
            idtA_shifted = idtA 
            timeLag = timedelta(seconds=self.catA['Lag_In_Track'][i])
            idtB_shifted = np.where(
                (self.tenHzB['dateTime'] >= t0-dt+timeLag-overlapW) &  
                (self.tenHzB['dateTime'] < t0+dt+timeLag+overlapW) &
                (self.tenHzB['dos1rate'] != -1E31)
                )[0]
            t_sA = t0
            t_sB = t0 + timeLag
        else:
            idtB_shifted = idtB 
            timeLag = timedelta(seconds=self.catB['Lag_In_Track'][i])
            idtA_shifted = np.where(
                (self.tenHzA['dateTime'] >= t0-dt-timeLag-overlapW) &  
                (self.tenHzA['dateTime'] < t0+dt-timeLag+overlapW) &
                (self.tenHzA['dos1rate'] != -1E31)
                )[0]
            t_sA = t0 - timeLag
            t_sB = t0
        return idtA, idtB, idtA_shifted, idtB_shifted, t0, t_sA, t_sB


    def _load_catalog(self, sc_id, catPath, v):
        """ 
        Loads the microburst catalog of version v spacecraft given 
        by self.sc_id. 
        """
        fPath = os.path.join(catPath, 
            'AC6{}_microbursts_v{}.txt'.format(sc_id.upper(), v))
        print('Loading AC6{}_microbursts_v{}.txt catalog'.format(
            sc_id.upper(), v))
        catalog = pd.read_csv(fPath, na_values=-1E31)
        # Convert the catalog times to datetime objects
        catalog['dateTime'] = pd.to_datetime(catalog['dateTime'])
        return catalog


if __name__ == '__main__':
    import time
    start_time = time.time()
    catalog_version = 4
    c = CumulativeDist(catalog_version, cc_width=2)
    try:
        c.loop()
    finally:
        saveDir = './../data/coincident_microbursts_catalogues'
        saveName = 'AC6_coincident_microbursts_2scc_v9.txt'
        c.save_catalog(savePath=os.path.join(saveDir, saveName))
        print("Ran in {} s".format(time.time() - start_time))
        #pass
