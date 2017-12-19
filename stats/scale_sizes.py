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

    def scaleSizeHist(self, bins, range=None, norm=None, visualizeNorm=True,
                      lowAE=None, highAE=None):
        """
        This function plots the histogram of microburst scale sizes.
        """
        # Filter out errororus separations
        validInd = np.where(self.dataDict['Dist_Total'] != -1E31)[0]
        if lowAE is not None: # If lower bound AE is set.
            lowAEInd = np.where(self.dataDict['AE'] >= lowAE)[0]
            validInd = np.array(sorted(list(set(validInd) & set(lowAEInd))))
        if highAE is not None: # If lower bound AE is set.
            highAEInd = np.where(self.dataDict['AE'] <= highAE)[0]
            validInd = np.array(sorted(list(set(validInd) & set(highAEInd))))

        hist, bin_edges = np.histogram(self.dataDict['Dist_Total'][validInd],
            bins=bins)
        width = 0.9 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        if norm is not None:
            hist = np.divide(hist, norm)
        if visualizeNorm:
            fig, (ax, bx) = plt.subplots(2, sharex=True)
            # histNorm, bin_edges = np.histogram(self.count_norm,
            #    bins=self.dist_norm)
            #bx.bar(center, hist, align='center', width=width)
            bx.plot(self.dist_norm+width/2, self.count_norm/86400)
            bx.set(xlabel='Separation (km)', ylabel = 'Days at separation')
        else:
            fig, ax = plt.subplots(1)

        ax.bar(center, hist*86400, align='center', width=width)
        ax.set(ylabel='flashes/day', 
            title='AC6 {} scale size distribution'.format(self.burstType), 
            xlabel='Separation (km)')
        if visualizeNorm:
            return (ax, bx)
        else:
            return (ax, None)

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
            #next(reader) # Skip first line of header
            self.keys = next(reader)
            # Remove comment and empty space chars
            self.keys[0] = self.keys[0][2:] 
            # Empty data file
            data = np.nan*np.ones((N, len(self.keys)), dtype=object)
        return data

if __name__ == '__main__':
    #fPath = os.path.abspath('./../data/curtain_catalogues/curtains_catalogue.txt')
    fPath = os.path.abspath('./../data/flash_catalogues/flashes_catalogue.txt')
    normPath = os.path.abspath('./dist_norm.txt')
    ss = ScaleSize(fPath, burstType='Flashes')
    #ss.loadNormalization(normPath)
    #ax, bx = ss.scaleSizeHist(ss.dist_norm, norm=ss.count_norm[:-1])
    #ax, bx = ss.scaleSizeHist(ss.dist_norm)

    # Load normalization
    normDist = np.load(os.path.abspath('./normalization/dist_vals.npy'))
    normCounts = np.load(os.path.abspath('./normalization/dist_norm.npy'))

    # Rebin normalization histogram to n km (every nth point)
    n = 10
    normDist = normDist[::n]
    normCounts = np.array([np.sum(normCounts[i*n:(i+1)*n] ) for i in range(len(normCounts)//n)])

    ss.dist_norm = normDist[:-1]
    ss.count_norm = normCounts

    ax, _ = ss.scaleSizeHist(normDist, visualizeNorm=True, norm=normCounts)
    ax.set_xlim(0, 300)
    plt.show()