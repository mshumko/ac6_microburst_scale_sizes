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

        self._loadData() # Load microburst data

        if self.iType.lower() == 'kp':
            self._loadKp() # Load kp index
        elif self.iType.lower() == 'ae':
            self._loadAe()
        return

    def appendIndex(self, timeKey='dateTimeA'):
        """
        This function will loop over the microburst data set and look for
        matching geomagnetic indicies.
        """
        self.matchingIndex = np.nan*np.ones(len(self.dataDict['lat']), dtype=int)

        # Loop over microburst times
        for (i, t) in enumerate(self.dataDict[timeKey]):
            idx = np.where((t >= self.indexTimes[:-1]) & 
                (t < self.indexTimes[1:]))[0]
            print(t, idx, self.indexTimes[idx])
            if len(idx) == 0:
                break
            self.matchingIndex[i] = self.index[idx]
        print(self.matchingIndex)
        return

    def saveData(self):
        """
        This file saves the microburst data to test.txt for now...
        """
        with open(self.dataPath, 'w', newline='') as f:
            writer = csv.writer(f)

            # Header
            #self.keys[0] = '# {}'.format(self.keys[0])
            writer.writerow(
                np.concatenate((self.keys, [self.iType.upper()]) 
                ))
            # Loop over lines and save file.
            for (i, line) in enumerate(self.data):
                writer.writerow(
                    np.concatenate( (line, [self.matchingIndex[i]]) )
                    )
        return

    def _loadData(self):
        """
        This function will load in the microburst catalogues
        from both spacecraft
        """
        self.data = self._inputDataSkelletons(self.dataPath)

        with open(self.dataPath, 'r') as f:
            reader = csv.reader(f)
            next(reader) # Skip header
            next(reader)

            for i, line in enumerate(reader): # Read in the data
                self.data[i] = line
            # Now format data array into a dictionary for each spacecraft
            self.dataDict = {}
            for i, key in enumerate(self.keys):
                self.dataDict[key] = self.data[:, i]
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
            #self.keys[0] = self.keys[0][2:] 
            # Empty data file
            data = np.nan*np.ones((N, len(self.keys)), dtype=object)
        return data

    def _loadKp(self):
        """
        This function will load in the Kp data
        """
        self.indexTimes = np.nan*np.ones((0, 1), dtype=object)
        self.index = np.nan*np.ones((0, 1), dtype=int)

        # Loop over years in the data and load in that index.
        dataYears = sorted(set([i.year for i in self.dataDict['dateTimeA']]))
        for year in [2014]:
            fName = '{}_{}.txt'.format(year, self.iType.lower())
            indexData = spacepy.datamodel.readJSONheadedASCII(
                os.path.join(self.indexDir, fName))
            self.indexTimes = np.append(self.indexTimes, indexData['dateTime'])
            self.index = np.append(self.index, indexData[self.iType.lower()])
        # Convert dateTimes
        self.indexTimes = np.array([dateutil.parser.parse(t) for t in self.indexTimes])
        return

    def _loadAe(self, aeType='AE', headerSize=15):
        """

        """
        self.indexTimes = np.array([], dtype=object)
        self.index = np.array([], dtype=int)

        dTypes = ['AE', 'AU', 'AL', 'A0']
        idx = 3 + dTypes.index(aeType)
        
        dataYears = sorted(set([i.year for i in self.dataDict['dateTimeA']]))
        for year in dataYears:
            fName = '{}_{}.txt'.format(year, self.iType.lower())
            with open(os.path.join(self.indexDir, fName), 'r') as f:
                reader = csv.reader(f, skipinitialspace=True, delimiter=' ')
                data = np.array(list(reader))

                # Skip the header and convert times
                self.indexTimes = np.append(self.indexTimes, 
                    [dateutil.parser.parse(' '.join(line[0:2])) for line in 
                    data[headerSize:]])
                # Save indicies.
                indicies = list(map(float, [line[idx] for line in data[headerSize:]]))
                self.index = np.append(self.index, indicies)
        return

if __name__ == '__main__':
    iType = 'ae'
    dataType = 'flashes'
    indexDir = '/home/mike/research/geomag_indicies/ae'

    if dataType == 'flashes':
        dataPath = ('/home/mike/research/ac6-microburst-scale-sizes/data/'
                    'flash_catalogues/flash_catalogue_v2.txt')
    else:
        dataPath = ('/home/mike/research/ac6-microburst-scale-sizes/data/'
                    'curtain_catalogues/curtain_catalogue.txt')

    appendObj = AppendGeoMagIdx(iType, dataPath, indexDir)
    appendObj.appendIndex()
    appendObj.saveData()
