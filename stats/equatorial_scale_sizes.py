# This script will map the scale sizes to the magnetic equator.
import os
import csv
import numpy as np
import itertools
import dateutil.parser

import IRBEM

Re = 6371 # km

class EquatorScaleSize(IRBEM.MagFields):
    def __init__(self, cPath, Bmodel='OPQ77', fltDict=None):
        """

        """
        IRBEM.MagFields.__init__(self, kext=Bmodel)

        # Load and filter the catalog data
        self.cData = self._load_data_(cPath)
        if fltDict is not None:
            self._filter_data(fltDict)
        return

    def loop(self):
        """
        This method loops over the filtered catalog microbursts and calls 
        the mapToEquator function.
        """
        self.dEq = np.nan*np.ones(len(self.cData['burstType']))
        z = zip(self.cData['dateTimeA'], self.cData['lat'], self.cData['lon'],
                self.cData['alt'], self.cData['Dist_In_Track'])
                

        for i, (t, lat, lon, alt, d) in enumerate(z):
            X = {'dateTime':t, 'x1':alt, 'x2':lat, 'x3':lon}
            self.dEq[i] = self.mapToEquator(X, d)
        return

    def mapToEquator(self, X, d, maginput=None):
        """
        This function calculates the equatorial scale size in latitude, given a 
        2x3 array of lat, lon, alt.
        """
        if not hasattr(X['x1'], '__len__'):
            #len(X['x1']) == 1: # If 1d array.
            dLat = deltaLat(d, X['x1'])
            X2 = X.copy() # Copy X dict to modify the lat values in both.
            # Add half of delta latitude from the separation.
            X['x2'] -= dLat/2  
            X2['x2'] += dLat/2
        Xeq1 = self.find_magequator(X, maginput)['XGEO']
        Xeq2 = self.find_magequator(X2, maginput)['XGEO']
        d = Re*np.linalg.norm(Xeq1-Xeq2)
        return d

    def _load_data_(self, fPath):
        """
        This method reads in a csv file, and saves it to a dictionary.
        """
        with open(fPath) as f:
            r = csv.reader(f)
            keys = next(r)
            rawData = np.array(list(r))
        data = {}
        for (i, key) in enumerate(keys):
            data[key] = rawData[:, i]

        # Now convert the times array(s) to datetime, and all others except 'burstType' to a float.
        for key in itertools.filterfalse(lambda x:'dateTime' in x or 'burstType' in x, data.keys()):
                data[key] = data[key].astype(float)
        for key in filter(lambda x: 'dateTime' in x, data.keys()):
            data[key] = np.array([dateutil.parser.parse(i) for i in data[key]])
        return data

    def _filter_data(self, fltDict):
        """
        This methdo will filter the catalog. Filter dict (fltDict)
        must contain a key that is in the catalog file, and values 
        to filter to. If only one value is passed, equality will be
        checked, if two values are passed, a range of values will
        be found. 
        """
        ind = range(len(self.cData['burstType']))
        for (key, value) in fltDict.items():
            # Check if user supplied a list of values (and not a string).
            if hasattr(value, '__len__') and not isinstance(value, str):
                # Filter by range
                idx = np.where((self.cData[key] > np.min(value)) & 
                                (self.cData[key] < np.max(value)))[0]
            else:
                # Filter by equality
                idx = np.where(self.cData[key] == value)[0]
            ind = np.intersect1d(ind, idx)
        
        #idx = np.where(self.cData[key] == value)[0]
        for key in self.cData.keys():
            self.cData[key] = self.cData[key][ind]
        return

def deltaLat(d, alt):
    """
    This function calculates the latitude angle from an arc length 
    (effectively straight line) from a separation d, at an altitude alt.
    d and alt must be in units of km.
    """
    dLat = 180/np.pi*d/(Re+alt)
    return dLat

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from datetime import datetime

    Lrange = range(3, 8)

    for (lL, uL) in zip(Lrange[:-1], Lrange[1:]):
        cDir = '/home/mike/research/ac6-microburst-scale-sizes/data/flash_catalogues'
        cName = 'flash_catalogue_v2_sorted.txt'
        ss = EquatorScaleSize(os.path.join(cDir, cName), 
                            fltDict={'burstType':'flash', 'Lm_OPQ':[lL, uL]})
        ss.loop()

        # Visualize these scale sizes
        import matplotlib.pyplot as plt
        plt.hist(ss.dEq, bins=np.arange(0, 2001, 100))
        plt.xlabel('Equatorial Scale Size [km]')
        plt.ylabel('Count')
        plt.title('Equatorial Microburst Scale Sizes | OPQ model | {} < L < {}'.format(lL, uL))
        plt.savefig('/home/mike/Dropbox/0_grad_work/ac6-flashes-curtains/plots'
                    '/{}/equatorial_microburst_scale_sizes_{}_L_{}.png'.format(
                    datetime.now().date(), lL, uL))

        plt.close()