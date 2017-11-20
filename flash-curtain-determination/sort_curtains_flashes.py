# This program will sort microbursts by curtains and flashes and output the
# data to the ./../data/<flash or curtain>_catalogues.
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import dateutil.parser
import itertools 
import time
import csv
import os

class SortMicrobursts:
    def __init__(self, indir, date=None, flashesSwitch=True):
        self.date = date
        self._loadMicroburstCatalogues(indir)
        self.flashesSwitch = flashesSwitch
        return

    def simpleFindMatches(self, thresh=0.2):
        """
        Loop over list A and look for microbursts in list B at
        times +/- thresh. If In_Track_Lag < thresh, give 
        warning or error out.
        """
        self.flashes = np.nan*np.ones((0, len(self.keys)+2))
        dT = timedelta(seconds=thresh)
        
        for (i, tA) in enumerate(self.dataA['dateTime']):
            matchIdt = np.where(
                (tA >= self.dataB['dateTime'] - dT) & 
                (tA <= self.dataB['dateTime'] + dT))
            if len(matchIdt[0]) != 1: 
                continue # If 0 or > 1 matches found, ignore.
                
            # Save the same keys as input data except save
            # each spacecraft's times of flashes and their
            # dos1rates.
            # Calc stats on separation and location
            stats = []
            for key in self.keys[2:]:
                stats.append(np.mean([self.dataA[key][i], 
                    self.dataB[key][matchIdt[0][0]] ]))
            dataLine = np.array(
                [tA, self.dataA['dos1rate'][i], 
                self.dataB['dateTime'][matchIdt[0][0]],
                self.dataB['dos1rate'][matchIdt[0][0]] ])
            dataLine = np.concatenate((dataLine, stats))

            # Now save to a large 2-D array
            self.flashes = np.row_stack((self.flashes, dataLine))
        print(self.flashes)
        return

    def saveData(self, ):

        return

    def _loadMicroburstCatalogues(self, fDir):
        """
        This function will load in the microburst catalogues
        from both spacecraft
        """
        
        for sc_id in ['A', 'B']:
            if self.date is None:
                fName = 'AC6{}_microbursts.txt'.format(
                    sc_id)
            else:
                fName = 'AC6{}_{}_microbursts.txt'.format(
                    sc_id, self.date.date())
            data = self._inputDataSkelletons(sc_id, os.path.join(fDir, fName))

            with open(os.path.join(fDir, fName)) as f:
                reader = csv.reader(f)
                next(reader) # Skip header
                next(reader)

                for i, line in enumerate(reader): # Read in the data
                    data[i] = line
                # Now format data array into a dictionary for each spacecraft
                dataDict = {}
                for i, key in enumerate(self.keys):
                    dataDict[key] = data[:, i]
                # Loop over all keys but datetimes and set the array types to be floats.
                for key in itertools.filterfalse(lambda x: x == 'dateTime', self.keys):
                    dataDict[key] = dataDict[key].astype(float)

                # Convert datetimes
                dataDict['dateTime'] = np.array([dateutil.parser.parse(i) for i in dataDict['dateTime']])
                # Assign data to correct attribute
                if sc_id.upper() == 'A':
                    self.dataA = dataDict
                else:
                    self.dataB = dataDict
        return

    def _inputDataSkelletons(self, sc_id, fPath):
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
    startTime = time.time()
    inDir = os.path.abspath('./../data/daily_microburst_catalogues/')
    date = datetime(2017, 1, 11)
    sorter = SortMicrobursts(inDir, date=date)
    sorter.simpleFindMatches()
    print('Run time: {}'.format(round(time.time()-startTime)))
