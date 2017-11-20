# These programs will analyze the flash and curtain scale sizes.
import matplotlib.pyplot as plt
import dateutil.parser 
import numpy as np
import itertools
import csv
import sys
import os

class ScaleSize:
    def __init__(self, fPath, burstType='flashes'):
        self._loadData(fPath)
        self.burstType = burstType
        return

    def scaleSizeHist(self, bins, range=None, norm=None, visualizeNorm=True):
        """
        This function plots the histogram of microburst scale sizes.
        """
        validInd = np.where(self.dataDict['Dist_Total'] != -1E31)[0]
        hist, bin_edges = np.histogram(self.dataDict['Dist_Total'][validInd],
            bins=bins)
        print(bin_edges)
        width = 0.9 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        if norm is not None:
            hist = np.divide(hist, norm)
        if visualizeNorm:
            fig, (ax, bx) = plt.subplots(2, sharex=True)
            #histNorm, bin_edges = np.histogram(self.count_norm,
            #    bins=self.dist_norm)
            #bx.bar(center, histNorm, align='center', width=width)
            bx.plot(self.dist_norm+width/2, self.count_norm)
            bx.set(xlabel='Separation (km)', ylabel = 'Days',
            title='AC6B mean daily separation with 10 Hz data')
        else:
            fig, ax = plt.subplot()
        ax.bar(center, hist, align='center', width=width)
        ax.set(ylabel='# of flashes', 
            title='AC6 {} scale size distribution'.format(self.burstType) )
        if visualizeNorm:
            return ax, bx
        else:
            return ax

    def _loadData(self, fPath):
        """
        This function will load in the microburst catalogues
        from both spacecraft
        """
        data = self._inputDataSkelletons(fPath)

        with open(fPath) as f:
            reader = csv.reader(f)
            next(reader) # Skip header
            next(reader)

            for i, line in enumerate(reader): # Read in the data
                data[i] = line
            # Now format data array into a dictionary for each spacecraft
            self.dataDict = {}
            for i, key in enumerate(self.keys):
                self.dataDict[key] = data[:, i]
            # Loop over all keys but datetimes and set the array types to be floats.
            for key in itertools.filterfalse(lambda x:'dateTime' in x, self.keys):
                self.dataDict[key] = self.dataDict[key].astype(float)

            # Convert datetimes
            for key in filter(lambda x: 'dateTime' in x, self.keys):
                self.dataDict[key] = np.array([dateutil.parser.parse(i) for i in self.dataDict[key]])
        return

    def loadNormalization(self, fPath):
        """
        This function will load in the normalization data 
        """
        with open(fPath, 'r') as f:
            reader = csv.reader(f)
            next(reader) # Skip header
            next(reader)
            norm = np.array(list(reader), dtype=int)
        self.dist_norm = norm[:, 0]
        self.count_norm = norm[:, 1]

    def _inputDataSkelletons(self, fPath):
        """
        Create empty dictionary and arrays for microburst catalogues. 
        """
        # Scan the file to get number of lines.
        with open(fPath, 'r') as f:
            N = sum(1 for row in f) - 2

        # Get data keys and data from file
        with open(fPath, 'r') as f:
            reader = csv.reader(f)
            next(reader) # Skip first line of header
            self.keys = next(reader)
            # Remove comment and empty space chars
            self.keys[0] = self.keys[0][2:] 
            # Empty data file
            data = np.nan*np.ones((N, len(self.keys)), dtype=object)
        return data

if __name__ == '__main__':
    fPath = os.path.abspath('./../data/curtain_catalogues/curtains_catalogue.txt')
    #fPath = os.path.abspath('./../data/flash_catalogues/flashes_catalogue.txt')
    normPath = os.path.abspath('./dist_norm.txt')
    ss = ScaleSize(fPath, burstType='Curtains')
    ss.loadNormalization(normPath)
    ax, bx = ss.scaleSizeHist(ss.dist_norm, norm=ss.count_norm[:-1])
    ax.set_xlim(0, 150)
    plt.show()