import itertools
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

from mission_tools.ac6 import read_ac_data

Re=6371 # km

class SignificantNoise:
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
        
if __name__ == '__main__':
    sc_id = 'A'
    CC_thresh = 0.8
    date = datetime(2016, 10, 14)
    signif = SignificantNoise(sc_id, date)
    signif.main()
    signif.calc_CDF(CC_thresh=CC_thresh)

    plt.pcolormesh(signif.XX, signif.YY, signif.cdf_grid)
    plt.title('AC6A random CC | {} | CC_thresh = {}'.format(date.date(), CC_thresh))
    plt.xlabel(r'$\Lambda$')
    plt.ylabel(r'$\Lambda$')
    plt.colorbar()
    plt.savefig('AC6_significant_random_CC.png', dpi=300)
    #plt.show()