import os
import numpy as np
import csv
import dateutil.parser
import itertools
import matplotlib.pyplot as plt

class GlobalDist:
    def __init__(self, cPath, normPath=None, normBinPath=None, fltDict=None):
        """
        This class bins the cataloged detections by L-MLT or L-lon. The 
        catalog can be filtered by the fltDict dictionary. 
        """
        self._load_cat_(cPath)
        if normPath is not None: # Load normalization
            self._load_norm_(normPath, normBinPath)
        
        if fltDict is not None:
            self._filter_data(fltDict) 
        return

    def bin_detections(self, keys=None):
        """
        This method will bin the detections from the catalog.
        """
        if keys is None:
            self.H, _, _ = np.histogram2d(self.cData[self.normKeys[0]], 
                                    self.cData[self.normKeys[1]], 
                                    bins=self.normBins)
        else:
            self.H, _, _ = np.histogram2d(self.cData[keys[0]], 
                                    self.cData[keys[1]], bins=self.normBins)
        if hasattr(self, 'norm'):
            self.H = np.nan_to_num(self.H/self.norm)
            #self.H /= self.norm
        return

    def make_L_MLT_map(self):
        f, ax = plt.subplots(2)
        im = ax[0].pcolormesh(self.normBins[1], self.normBins[0], self.H)
        ax[0].set(xlabel=self.normKeys[1], ylabel=self.normKeys[0], title='Microbururst Occurance Rates')
        plt.colorbar(im, ax=ax[0], label = 'microbursts/s')

        im = ax[1].pcolormesh(self.normBins[1], self.normBins[0], self.norm)
        ax[1].set(xlabel=self.normKeys[1], ylabel=self.normKeys[0], title='Normalization')
        plt.colorbar(im, ax=ax[1], label = 'Seconds/bin')

        plt.tight_layout()
        return


    def _load_cat_(self, fPath):
        """
        This method reads in a csv file, and saves it to a dictionary.
        """
        with open(fPath) as f:
            r = csv.reader(f)
            keys = next(r)
            rawData = np.array(list(r))
        self.cData = {}
        for (i, key) in enumerate(keys):
            self.cData[key] = rawData[:, i]

        # Now convert the times array(s) to datetime, and all others except 'burstType' to a float.
        for key in itertools.filterfalse(lambda x:'dateTime' in x or 'burstType' in x, self.cData.keys()):
                self.cData[key] = self.cData[key].astype(float)
        for key in filter(lambda x: 'dateTime' in x, self.cData.keys()):
            self.cData[key] = np.array([dateutil.parser.parse(i) for i in self.cData[key]])
        return

    def _load_norm_(self, normPath, binPath):
        """
        This method loads in the normalization and bin files for the global maps
        """
        with open(normPath) as f:
            r = csv.reader(f)
            self.normKeys = next(r)
            self.norm = np.array(list(r), dtype=float)

        with open(binPath) as f:
            r = csv.reader(f)
            next(r) # Skip header.
            self.normBins = [[float(x) for x in y] for y in list(r)]
        return


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

if __name__ == '__main__':
    cDir = '/home/mike/research/ac6-microburst-scale-sizes/data/flash_catalogues'
    cName = 'flash_catalogue_v2_sorted.txt'
    # Normalization dir and name
    normDir = '/home/mike/research/ac6-microburst-scale-sizes/data/norm'
    normName = 'ac6_L_lon_norm.csv'
    normBinName = 'ac6_L_lon_bins.csv'
    
    ss = GlobalDist(os.path.join(cDir, cName), os.path.join(normDir, normName), 
                    os.path.join(normDir, normBinName), 
                    fltDict={'burstType':'flash'})
    ss.bin_detections()
    ss.make_L_MLT_map()
    plt.show()
