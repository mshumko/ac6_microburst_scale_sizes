import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

from mission_tools.ac6 import read_ac_data

Re=6371 # km

class SignificantNoise:
    def __init__(self, sc_id, date, microburst_time, CC_width=1):
        """
        This class calculates the fraction of random cross-correlations 
        above CC_threshold. This gives an indication of the fraction of
        signigicant detections due to only random noise.
        """
        self.sc_id = sc_id
        self.date = date
        self.CC_width = CC_width
        self.microburst_time = microburst_time
        self._load_data() # Load 10 Hz data from date.
        return

    def main(self, N=1E6, bins=np.arange(0, 70, 5)):
        """
        This method loops N times and for each iteration picks a random
        CC_width dos1 count rates and cross correlates it against the 
        CC_width count rates centered on self.micrburst_time
        """
        #self._find_reference_microbirst_counts()
        #self._calc_dist_relative_to_ref()

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

    # def _calc_dist_relative_to_ref(self):
    #     """ 
    #     Calculate the distance from the position at self.microburst_time 
    #     to all other data points.
    #     """
    #     idt = np.where(self.tenHzData['dateTime'] > self.microburst_time)[0][0]
    #     ref_pos = np.array([self.tenHzData['lat'][idt],
    #                         self.tenHzData['lon'][idt],
    #                         self.tenHzData['alt'][idt]])
    #     X1 = np.tile(ref_pos, (len(self.tenHzData['lat']), 1))
    #     X2 = np.array([self.tenHzData['lat'], 
    #                     self.tenHzData['lon'],
    #                     self.tenHzData['alt']]).T
    #     self.dist = self.haversine(X1, X2)
    #     return

    # def haversine(self, X1, X2):
    #     """
    #     Implementation of the haversine foruma to calculate total distance
    #     at an average altitude. X1 and X2 must be N*3 array of 
    #     lat, lon, alt.
    #     """
    #     X1 = np.asarray(X1)
    #     X2 = np.asarray(X2)
    #     R = (Re+(X1[:, 2]+X2[:, 2])/2)
    #     s = 2*np.arcsin( np.sqrt( np.sin(np.deg2rad(X1[:, 0]-X2[:, 0])/2)**2 + \
    #                     np.cos(np.deg2rad(X1[:, 0]))*np.cos(np.deg2rad(X2[:, 0]))*\
    #                     np.sin(np.deg2rad(X1[:, 1]-X2[:, 1])/2)**2 ))
    #     return R*s
    
    def _find_reference_microbirst_counts(self):
        """ 
        Find and return the dos1 count rates that are centered on 
        self.microburst_time.
        """
        tRange = [self.microburst_time - timedelta(seconds=self.CC_width/2),
                       self.microburst_time + timedelta(seconds=self.CC_width/2)]
        idt = np.where((self.tenHzData['dateTime'] > tRange[0]) & 
                       (self.tenHzData['dateTime'] <= tRange[1]))[0]
        self.ref_counts = self.tenHzData['dos1rate'][idt]
        return 

    def _load_data(self):

        self.tenHzData = read_ac_data.read_ac_data_wrapper(
                            self.sc_id, self.date)
        return
        
if __name__ == '__main__':
    sc_id = 'A'
    date = datetime(2016, 10, 14)
    microburst_time = datetime(2016, 10, 14, 4, 27, 13)
    signif = SignificantNoise(sc_id, date, microburst_time)
    signif.main()
