# This file appends Kp index to the curtain and flash catalogues.
import sys
import os
import csv
import numpy as np
import glob
import dateutil.parser
import itertools

import spacepy.datamodel

class AppendGeoMagIdx:
    """

    """
    def __init__(self, iType, dataPath, indexDir):
        """
        This class will read in the microburst catalogues, 
        and append a geomagnetic index to the end of it.
        """
        self.iType = iType
        self.indexDir = indexDir
        self.dataPath = dataPath

        # Load data and Kp
        self._loadData()
        self._loadIndicies()
        return

    def appendKp(self):
        
        return

    def _loadData(self):
        """
        This function will load in the microburst catalogues
        from both spacecraft
        """
        data = self._inputDataSkelletons(self.dataPath)

        with open(self.dataPath, 'r') as f:
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
                self.dataDict[key] = np.array([dateutil.parser.parse(i) 
                    for i in self.dataDict[key]])
        return

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

    def _loadIndicies(self):
        """
        This function will load in the Kp data
        """
        self.indexTimes = np.nan*np.ones((0, 1), dtype=object)
        self.index = np.nan*np.ones((0, 1), dtype=int)

        # Loop over years in the data and load in that index.
        dataYears = sorted(set([i.year for i in self.dataDict['dateTimeA']]))
        for year in dataYears:
            fName = '{}_{}.txt'.format(year, self.iType.lower())
            indexData = spacepy.datamodel.readJSONheadedASCII(
                os.path.join(self.indexDir, fName))
            self.indexTimes = np.append(self.indexTimes, indexData['dateTime'])
            self.index = np.append(self.index, indexData[self.iType.lower()])
        # Convert dateTimes
        self.indexTimes = np.array([dateutil.parser.parse(t) for t in self.indexTimes])
        return

if __name__ == '__main__':
    iType = 'Kp'
    dataType = 'flashes'
    indexDir = '/home/mike/research/firebird/data_processing/geomag_indicies/indicies'

    if dataType == 'flashes':
        dataPath = ('/home/mike/research/ac6-microburst-scale-sizes/data/'
                    'flash_catalogues/flashes_catalogue.txt')
    else:
        dataPath = ('/home/mike/research/ac6-microburst-scale-sizes/data/'
                    'curtain_catalogues/curtains_catalogue.txt')

    appendObj = AppendGeoMagIdx(iType, dataPath, indexDir)