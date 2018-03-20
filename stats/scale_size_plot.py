# This code will generate the scale size distribution histograms
import matplotlib.pyplot as plt
import dateutil.parser 
import numpy as np
import itertools
import csv
import sys
import os

class ScaleSizeDist:
    def __init__(self, cPath, nPath, fltDict=None):
        """
        Load the catalog file located at "cPath" and save it to a cData 
        attribute of the class, and load the normalization file from
        nPath, and save it to a nData attribute.
        """
        self.cData = self._load_data_(cPath)
        self.nData = self._load_data_(nPath)
        
        if fltDict is not None:
            self._filter_data(fltDict) 
        return

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

    def plot_hist(self, ax=None, norm=True):
        if ax is None:
            fig, self.ax = plt.subplots()
        else:
            self.ax = ax

        bins = self.nData['Separation [km]']
        #bins = np.append(bins, (2*bins[-1] - bins[-2]))
        if norm: # Normalize the histogram
            H, bin_edges = np.histogram(self.cData['Dist_Total'],
                bins=np.append(bins, (2*bins[-1] - bins[-2])))
            width = 0.9 * (bins[1] - bins[0])
            center = (bins[:-1] + bins[1:]) / 2
            H = 86400*np.divide(H, self.nData['Seconds'])
            self.ax.bar(bins+width/2, H, align='center', width=width)
            self.ax.set(ylabel='flash/day')
        else:
            self.ax.hist(self.cData['Dist_Total'], bins=bins)

        if ax is None:
            plt.show()
        return self.ax

    def plot_norm(self, bx=None):
        """
        This method will plot the normalization coefficients
        """
        if bx is None:
            fig, self.bx = plt.subplots()
        else:
            self.bx = bx
        bins = self.nData['Separation [km]']
        width = 0.9 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        self.bx.bar(bins+width/2, self.nData['Seconds']/86400, align='center', width=width)
        self.bx.set(xlabel='Separation (km)', ylabel = 'days')
        return self.bx

if __name__ == '__main__':
    # Catalog dir and name
    cDir = '/home/mike/research/ac6-microburst-scale-sizes/data/flash_catalogues'
    cName = 'flash_catalogue_v2_sorted.txt'
    # Normalization dir and name
    nDir = '/home/mike/research/ac6-microburst-scale-sizes/data/norm'
    nName = 'ac6_norm_all.txt'
    
    ss = ScaleSizeDist(os.path.join(cDir, cName), os.path.join(nDir, nName), 
                    fltDict={'burstType':'flash'})

    # Plot data
    fig, ax = plt.subplots(2, sharex=True)
    ss.plot_hist(ax=ax[0])
    ax[0].set_xlim(0, 100)
    ss.plot_norm(bx=ax[1])
    ax[0].set_title('Flash Scale Sizes')
    ax[1].set_title('Normalization')

    plt.savefig(('/home/mike/Dropbox/0_grad_work/ac6-flashes-curtains/'
                'plots/2018-03-14_scale_sizes/'
                'flash_scale_size_all.png'))
