import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
import dateutil.parser
from datetime import datetime, timedelta
import os
import csv
import glob

from mission_tools.ac6 import read_ac_data
from mission_tools.misc import obrienBaseline

Re=6371 # km

class SignificantNoise:
    """
    First iteration at solving this problem. Here I CCd random times and 
    binned them up by Lambda (invariant latitude) for each day. This shows
    an occurance rate of being too low as I would expect. This may be due 
    to CC random times and not microburst detections and not accounting for
    CCs at times with similar L-MLT-AE. 
    """
    def __init__(self, sc_id, date, CC_width=1, bins=np.arange(60, 70, 0.5)):
        """
        This class calculates the fraction of random cross-correlations 
        above CC_threshold. This gives an indication of the fraction of
        signigicant detections due to only random noise.
        """
        self.sc_id = sc_id
        self.date = date
        self.CC_width = CC_width
        self.bins = bins
        self._load_data() # Load 10 Hz data from date.
        return

    def main(self, N=50000, CC_width_thresh=1, verbose=True):
        """
        This method loops over a meshgrid specified by self.bins and
        for each bin calculates N CC's between random dos1 time series
        of length self.CC_width. CC_width_thresh is the wiggle room to 
        calculate the max CC.
        """
        self.XX, self.YY = np.meshgrid(self.bins, self.bins)
        self.CC_arr = np.nan*np.zeros((*self.XX.shape, N), dtype=float)
        lx, ly = self.XX.shape
        CC_data_width = 5/self.CC_width 
        InvLat_OPQ = np.abs(self.tenHzData['InvLat_OPQ'])

        # Loop over the meshgrid.
        for (i, j) in itertools.product(range(lx-1), range(ly-1)):
            bin_i_edges = (self.XX[i, j], self.XX[i, j+1])
            bin_j_edges = (self.YY[i, j], self.YY[i+1, j])

            # Find all data indicies in  InvLat bin
            idx = np.where((InvLat_OPQ > bin_i_edges[0]) & 
                            (InvLat_OPQ <= bin_i_edges[1]))[0]
            jdx = np.where((InvLat_OPQ > bin_j_edges[0]) & 
                            (InvLat_OPQ <= bin_j_edges[1]))[0]
            if verbose:
                print('Invlat Edges = ', bin_i_edges, bin_j_edges)
                print('Number of data points = ', len(idx), len(jdx))

            # Pick N random center CC indicies
            centers_i = np.random.choice(idx, size=N)
            centers_j = np.random.choice(jdx, size=N)

            for k in range(N):
                cc_ind_i = np.arange(centers_i[k] - CC_data_width, 
                                     centers_i[k] + CC_data_width, 
                                     dtype=int)
                cc_ind_j = np.arange(centers_j[k] - CC_data_width - CC_width_thresh, 
                                     centers_j[k] + CC_data_width + CC_width_thresh,
                                     dtype=int)
                if ((max(cc_ind_i) > len(InvLat_OPQ)-1 or max(cc_ind_j) > len(InvLat_OPQ)-1) or 
                    (min(cc_ind_i) < 0 or min(cc_ind_j) < 0)):
                    continue
                self.CC_arr[i, j, k] = self.CC(cc_ind_i, cc_ind_j)
        return

    def calc_CDF(self, CC_thresh=0.8):
        """
        For each bin on self.CC_arr, calculate the fraction of CC values greater
        than CC_thresh
        """
        self.cdf_grid = np.nan*np.zeros((len(self.bins)-1, len(self.bins)-1), dtype=float)

        lx, ly, _ = self.CC_arr.shape
        #idx = np.where(signif.CC_arr > 0.8)
        # Loop over the CC_arr meshgrid.
        for (i, j) in itertools.product(range(lx-1), range(ly-1)):
            ids = np.where(self.CC_arr[i, j, :] > CC_thresh)[0]
            idtotal = np.where(~np.isnan(self.CC_arr[i, j, :]))[0]
            self.cdf_grid[i, j] = len(ids)/len(idtotal)
        return

    def CC(self, iA, iB):
        """ 
        This method calculates the normalized cross-correlation 
        between two AC6 time series indexed by iA and iB.
        """
        norm = np.sqrt(len(self.tenHzData['dos1rate'][iA])*\
                       len(self.tenHzData['dos1rate'][iB])*\
                       np.var(self.tenHzData['dos1rate'][iA])*\
                       np.var(self.tenHzData['dos1rate'][iB])) 
        # Mean subtraction.
        x = (self.tenHzData['dos1rate'][iA] - 
            self.tenHzData['dos1rate'][iA].mean() )
        y = (self.tenHzData['dos1rate'][iB] - 
            self.tenHzData['dos1rate'][iB].mean() )
        # Cross-correlate
        ccArr = np.correlate(x, y, mode='valid')
        # Normalization
        ccArr /= norm
        return max(ccArr)

    def _load_data(self):

        self.tenHzData = read_ac_data.read_ac_data_wrapper(
                            self.sc_id, self.date)
        return

class SignificantFraction:
    """
    Second iteration at solving this problem. Here I CCd 10 random microbursts
    picked from each day against random time series in 4 < L < 8 on that day.
    No binning was attempted, and this class generates a CCtFrac array that
    contains the fraction of random events above a significance threshold
    as a function of day. This still shows a random chance rate that is
    about 5 times lower than I'd expect.
    """
    def __init__(self, sc_id, catalog_version, CC_window=10):
        """

        """
        self.sc_id = sc_id
        self.CC_window = CC_window
        catPath = ('/home/mike/research/ac6_microburst_scale_size'
                    's/data/microburst_catalogues')
        self._load_catalog(catPath, catalog_version)
        return

    def main(self, CC_window_thresh=1, CC_thresh=0.8, fPath=None, 
             N_microburst=10, N_CC=100):
        """ 
        Main method to loop over AC6 data and cross-correlate random
        time series.
        """
        # if fPath is None: 
        #     fPath = '/home/mike/research/ac6/ac6{}/ascii/level2'.format(self.sc_id)
        # paths = glob.glob(os.path.join(fPath, '*10Hz*.csv'))
        
        self._unique_dates()
        self.CCtFrac = np.nan*np.zeros(len(self.cat_dates), dtype=float)

        for i, date in enumerate(self.cat_dates):
            self.dayCC = np.nan*np.zeros((N_microburst, N_CC))
            # Load AC6 10 Hz data
            self._load_10Hz_data(date)
            # Pick N_microbursts number of microbursts (if there are that many.)
            self._find_microburst_times(date, N_microburst)
            
            # To prepare the CC, find dos1rate indicies that are in the 
            # rad belts.
            self.tenHzData['Lm_OPQ'] = np.abs(self.tenHzData['Lm_OPQ'])
            iBelt = np.where( (self.tenHzData['Lm_OPQ'] > 5) & 
                              (self.tenHzData['Lm_OPQ'] < 5.5) &
                              (self.tenHzData['flag'] == 0) )[0]
            if not len(iBelt):
                continue
            # Get N_CC rad belt indicies as centers for CC windows.
            jRandom = np.random.choice(iBelt, size=N_CC)

            # For each microburst, CC N_CC times.
            for row, iRow in enumerate(self.microburst_idx):
                for col, iCol in enumerate(jRandom):
                    iA, iB = self._get_CC_indicies(iRow, iCol, CC_window_thresh)
                    self.dayCC[row, col] = self.CC(iA, iB)
            # If dayCC array has at least one non-nan value.
            numerator = len(np.where(self.dayCC > CC_thresh)[0])
            denominator = len(np.where(~np.isnan(self.dayCC))[0])
            if denominator:
                self.CCtFrac[i] = (numerator/denominator)
        return

    def CC(self, iA, iB):
        """ 
        This method calculates the normalized cross-correlation 
        between two AC6 time series indexed by iA and iB.
        """
        norm = np.sqrt(len(self.tenHzData['dos1rate'][iA])*\
                       len(self.tenHzData['dos1rate'][iB])*\
                       np.var(self.tenHzData['dos1rate'][iA])*\
                       np.var(self.tenHzData['dos1rate'][iB])) 
        # Mean subtraction.
        x = (self.tenHzData['dos1rate'][iA] - 
            self.tenHzData['dos1rate'][iA].mean() )
        y = (self.tenHzData['dos1rate'][iB] - 
            self.tenHzData['dos1rate'][iB].mean() )
        # Cross-correlate
        ccArr = np.correlate(x, y, mode='valid')
        # Normalization
        ccArr /= norm
        return max(ccArr)
        
    def _get_CC_indicies(self, icA, icB, CC_window_thresh):
        """ 
        For center indicies iA and iB, return the CC indicies with a 
        CC_window_thresh threshold (in data points).
        """
        iAmin = icA-self.CC_window//2
        iAmax = icA+self.CC_window//2
        
        iBmin = icB-self.CC_window//2-CC_window_thresh
        iBmax = icB+self.CC_window//2+CC_window_thresh
        
        # Now check and fix cases where we are out of bounds.
        if iAmin < 0: iAmin = 0
        if iBmin < 0: iBmin = 0
        if iAmax >= len(self.tenHzData['dos1rate']): 
            iAmax = len(self.tenHzData['dos1rate'])-1
        if iBmax >= len(self.tenHzData['dos1rate']): 
            iBmax = len(self.tenHzData['dos1rate'])-1
           
        iA = np.arange(iAmin, iAmax)
        iB = np.arange(iBmin, iBmax)
        return iA, iB

    def _unique_dates(self):
        allDates_num = date2num(self.cat['dateTime']).astype(int)
        unique_dates_num = set(allDates_num)
        sorted_dates_num = sorted(list(unique_dates_num))
        self.cat_dates = num2date(sorted_dates_num)
        return 

    def _load_10Hz_data(self, date):
        """ Wrapper to load AC6 10 Hz data """
        print('Loading AC6-{} data from {}'.format(
            self.sc_id.upper(), date.date()))
        self.tenHzData = read_ac_data.read_ac_data_wrapper(
                            self.sc_id, date)
        return

    def _find_microburst_times(self, date, N):
        num_date = date2num(date)
        num_cat_dates = date2num(self.cat['dateTime']).astype(int)

        idM = np.where(num_date == num_cat_dates)[0]
        #print('idM', idM)
        # In case there are less than N total microbursts in the 
        # catalog on that day.
        if len(idM) < N: N = len(idM)
        # Choose N microburst indicies without replacement.
        self.microburst_idx = np.random.choice(idM, N, replace=False)
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

class CrossCorrelateMicrobursts(SignificantFraction):
    """
    This is the third iteration of the significance chance
    coincidence code. Here I loop over all of the detections
    from the catalog and save a time series snippet containing
    the microburst as well as time, L, MLT, and AE to a file. 
    This file is then read back in and then each micrburst
    is CCd against other microbursts observed in a similar
    L-MLT-AE bin.
    """
    def __init__(self, sc_id, catalog_version, CC_window=10, 
                timeSeriesPath='microburst_counts.csv', CC_thresh=0.8):
        self.timeSeriesPath = timeSeriesPath
        self.CC_thresh = CC_thresh
        super().__init__(sc_id, catalog_version, CC_window=10)
        return

    def saveTimeSeries(self, CC_window_thresh=1):
        """ 
        Main method to loop over AC6 data and save each microburst's
        the timeseries to a csv file.
        """
        self._unique_dates()

        if not os.path.isfile(self.timeSeriesPath):
            with open(self.timeSeriesPath, 'w') as f:
                w = csv.writer(f)
                w.writerow(['dateTime', 'Lm_OPQ', 'MLT_OPQ', 'AE']+
                        (self.CC_window + 2*CC_window_thresh)*['CountsArray'])

        for i, date in enumerate(self.cat_dates):
            # Load AC6 10 Hz data
            self._load_10Hz_data(date)
            # Find the baseline count rates.
            #validCounts = np.where(self.tenHzData['dos1rate'] > 0)[0]
            #self.baseline = obrienBaseline.obrienBaseline(
            #    self.tenHzData['dos1rate'][validCounts], cadence=0.1) 
            # Find the microbursts on this day
            self._get_daily_microbursts(date)
            with open(self.timeSeriesPath, 'a') as f:
                w = csv.writer(f)
                for j in self.daily_microbursts:
                    # Write the count rates to file here.

                    tenHzIdx = np.where(
                        self.cat['dateTime'][j] == self.tenHzData['dateTime']
                                        )[0]
                    assert len(tenHzIdx) == 1, ('Zero or > 1 '
                                                'microburst matches '
                                                'found! tenHzIdx={}'.format(
                                                    tenHzIdx))
                    # dos1 rate indicies to save
                    idc_start = tenHzIdx[0] - self.CC_window//2 - CC_window_thresh
                    idc_end   = tenHzIdx[0] + self.CC_window//2 + CC_window_thresh
                    # Skip the microburst if the count data is at the dos1rate 
                    # file boundary
                    if idc_start < 0 or idc_end >= len(self.tenHzData['dateTime']):
                        continue

                    # Skip the microburst if there are error values in the 
                    # count data.
                    if np.min(self.tenHzData['dos1rate'][idc_start:idc_end]) < 0:
                        continue

                    # Write data to a file.
                    w.writerow(
                        np.concatenate((
                            [date2num(self.cat['dateTime'][j])],
                            [self.cat['Lm_OPQ'][j]], 
                            [self.cat['MLT_OPQ'][j]],
                            [self.cat['AE'][j]],
                            self.tenHzData['dos1rate'][idc_start:idc_end]
                        ))
                    )
        return

    def binMicrobursts(self, L_bins=np.arange(4, 8, 0.5), 
                            MLT_bins=np.arange(0, 15, 0.5), 
                            AE_bins=np.arange(0, 3000, 100),
                            N_CC = 100):
        # Create 3d meshgrid to loop over
        LL, MLTMLT, AEAE = np.meshgrid(L_bins, MLT_bins, AE_bins)
        F = np.nan*np.zeros_like(LL)

        lL, lMLT, lAE = LL.shape

        # A nested loop that is three for loops deep.
        for i, j, k in itertools.product(range(lL), range(lMLT), range(lAE)):
            # Search for all microbursts in that bin. The 3d array indicies
            # steps were found through trial and error.
            iBursts = np.where(
                (self.d[:, 1] > LL[i, j, k]) & (self.d[:, 1] < LL[i, j+1, k]) &
                (self.d[:, 2] > MLTMLT[i+1, j, k]) & (self.d[:, 2] < MLTMLT[i+1, j, k]) &
                (self.d[:, 3] > AEAE[i, j, k]) & (self.d[:, 3] < AEAE[i, j, k+1])
            )[0]

            # Skip bins with only a "few" microbursts in this bin.
            if len(iBursts) < N_CC:
                continue

            # Pick N_CC random microbursts from the iBursts list to cross-correlate against.
            iBursts2 = np.random.choice(iBursts, size=N_CC)
            # Loop over the microbursts in that bin and cross-correlate them.
            CCarr = np.nan*np.zeros((len(iBursts), N_CC))

            for a, iBurst in enumerate(iBursts):
                # Now cross correlate iBurst agaianist the random microbursts in iCC
                for b, iBurst2 in enumerate(iBursts2):
                    CCarr[a, b] = self.CC(self.d[iBurst, 4:], self.d[iBurst2, 4:])
            
            # Now calculate the ratio of CCs above the threshold against all other CCs.
            numerator = len(np.where(CCarr > self.CC_thresh)[0])
            denominator = len(np.where(~np.isnan(CCarr))[0])
            # Avoid division by 0
            if denominator:
                F[i, j, k] = numerator/denominator
        return


    def CC(self, cA, cB):
        """ 
        This method calculates the normalized cross-correlation 
        between two AC6 time series indexed by iA and iB.
        """
        norm = np.sqrt(len(cA)*len(cB)*np.var(cA)*np.var(cB))
        # Mean subtraction.
        x = (cA - cA.mean())
        y = (cB - cB.mean())
        # Cross-correlate
        ccArr = np.correlate(x, y, mode='valid')
        # Normalization
        ccArr /= norm
        return max(ccArr)

    def _load_count_data(self):
        """ Load the microburst data saved to self.timeSeriesPath """
        self.d = np.genfromtxt(self.timeSeriesPath, delimiter=',', skip_header=1)
        return


    def _get_daily_microbursts(self, date):
        """
        Get microburst catalog indicies from date.
        """
        num_date = date2num(date)
        num_cat_dates = date2num(self.cat['dateTime']).astype(int)
        self.daily_microbursts = np.where(num_date == num_cat_dates)[0]
        return

if __name__ == '__main__':
    # sf = SignificantFraction('a', 5)
    # sf.main()
    # np.savetxt('random_CC_significance.csv', sf.CCtFrac)

    # # Visualize histogram of the cross-correlations.
    # frac = sf.CCtFrac[np.where(~np.isnan(sf.CCtFrac))[0]]
    # plt.hist(100*frac)
    # plt.show()

    ccm = CrossCorrelateMicrobursts('a', 5)
    ccm.saveTimeSeries()


# if __name__ == '__main__':
#     sc_id = 'A'
#     CC_thresh = 0.8
#     date = datetime(2015, 4, 25)
#     signif = SignificantNoise(sc_id, date)
#     signif.main()
#     signif.calc_CDF(CC_thresh=CC_thresh)

#     plt.pcolormesh(signif.XX, signif.YY, signif.cdf_grid)
#     plt.title('AC6A random CC | {} | CC_thresh = {}'.format(date.date(), CC_thresh))
#     plt.xlabel(r'$\Lambda$')
#     plt.ylabel(r'$\Lambda$')
#     plt.colorbar()
#     plt.savefig('AC6_significant_random_CC.png', dpi=300)
#     #plt.show()
