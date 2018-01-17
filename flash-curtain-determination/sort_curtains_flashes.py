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
    def __init__(self, inDir, outDir, date=None, flashesFlag=True):
        self.date = date
        self.inDir = inDir
        self.flashesFlag = flashesFlag
        self.outDir = outDir
        self._loadMicroburstCatalogues()
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

    def saveData(self, outName=None):
        """
        This function will save the curtain or flashes database to a csv file.
        """
        if self.flashesFlag:
            label = 'flash'
        else:
            label = 'curtain'

        if (outName is None) and (self.date is None):
            outName = '{}_catalogue.txt'.format(label)
        elif (outName is None) and (self.date is not None): 
            outName = '{}_{}_catalogue.txt'.format(self.date.date(), label)

        outKeys = np.concatenate((['dateTimeA', 'dos1rateA', 'dateTimeB', 'dos1rateB'], 
            self.keys[2:]))

        with open(os.path.join(self.outDir, outName), 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['{} catalogue created on {}'.format(
                label, datetime.now().date())])
            writer.writerow(outKeys)

            for line in self.flashes:
                writer.writerow(line)
        return

    def _loadMicroburstCatalogues(self):
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
            data = self._inputDataSkelletons(sc_id, os.path.join(self.inDir, fName))

            with open(os.path.join(self.inDir, fName)) as f:
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

        # Now shift timestamps by the in track lag for AC6-A if self.flashesFlag is False
        if not self.flashesFlag:
            self._shiftTimes()
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
            #next(reader) # Skip first line of header
            self.keys = next(reader)
            # Remove comment and empty space chars
            #self.keys[0] = self.keys[0][2:] 
            # Empty data file
            data = np.nan*np.ones((N, len(self.keys)), dtype=object)
        return data

    def _shiftTimes(self):
        """
        This internal method will shift the datetimes of one spacecraft
        by the Lag_In_Track
        """
        validLags = np.where(self.dataA['Lag_In_Track'] != -1E31)[0]
        self.dataA['dateTime'] = np.array([dT + timedelta(seconds=dS) for 
                                (dT, dS) in zip(self.dataA['dateTime'][validLags], 
                                self.dataA['Lag_In_Track'][validLags])])
        return

if __name__ == '__main__':
    burstType = 'curtain'
    startTime = time.time()
    inDir = os.path.abspath('./../data/microburst_catalogues/')
    outDir = os.path.abspath('./../data/{}_catalogues/'.format(burstType))

    if 'flash' in burstType:
        flashesFlag = True
    else:
        flashesFlag = False
    sorter = SortMicrobursts(inDir, outDir, flashesFlag=flashesFlag)
    sorter.simpleFindMatches()
    sorter.saveData()
    print('Run time: {}'.format(round(time.time()-startTime)))