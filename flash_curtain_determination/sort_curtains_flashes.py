# This program will sort microbursts by curtains and flashes and output the
# data to the ./../data/<flash or curtain>_catalogues.
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import dateutil.parser
import scipy.signal
import itertools 
import time
import csv
import sys
import os

sys.path.append('/home/mike/research/mission-tools/ac6/')
import read_ac_data

class SortMicrobursts:
    def __init__(self, inDir, outDir, date=None, flashesFlag=True, cadence=0.1):
        self.date = date
        self.inDir = inDir
        self.flashesFlag = flashesFlag
        self.outDir = outDir
        self._loadMicroburstCatalogues()
        self.cadence = cadence
        return

    def simpleFindMatches(self, T_THRESH=0.2, CC_WITH=0.5):
        """
        Loop over list A and look for microbursts in list B at
        times +/- T_THRESH. If In_Track_Lag < T_THRESH, give 
        warning or error out.
        """
        self.flashes = np.nan*np.ones((0, len(self.keys)+3))
        dT = timedelta(seconds=T_THRESH)

        currentDate = datetime.min
        for (i, tA) in enumerate(self.dataA['dateTime']):
            matchIdt = np.where(
                (tA >= self.dataB['dateTime'] - dT) & 
                (tA <= self.dataB['dateTime'] + dT))
            if len(matchIdt[0]) != 1: 
                continue # If 0 or > 1 matches found, ignore.

            # My attempt at functional programming.
            # Whenever current date changes, load in the new data.
            if currentDate.date() != tA.date():
                print('Loading data from {}'.format(tA.date()))
                # Read in data to caculate the cross correlation
                dataA = read_ac_data.read_ac_data_wrapper('A', tA)
                dataB = read_ac_data.read_ac_data_wrapper('B', tA)
                currentDate = tA

            # Calculate the cross correlation coefficent.
            try:
                cc = self._calcCrossCorr(dataA, dataB, tA, 
                    self.dataB['dateTime'][matchIdt[0]], 
                    width=CC_WITH)  
            except ValueError as err:
                if str(err) == 'math domain error':
                    print(err)
                    continue
                else:
                    raise

            # Save the same keys as input data except save
            # each spacecraft's times of flashes and their
            # dos1rates.
            # Calc stats on separation and location
            stats = []
            for key in self.keys[2:]:
                stats.append(np.mean([self.dataA[key][i], 
                    self.dataB[key][matchIdt[0][0]] ]))
            stats.append(round(cc, 2)) # Append the cross correlation coefficient
            dataLine = np.array(
                [tA, self.dataA['dos1rate'][i], 
                self.dataB['dateTime'][matchIdt[0][0]],
                self.dataB['dos1rate'][matchIdt[0][0]] ])
            dataLine = np.concatenate((dataLine, stats))

            # Now save to a large 2-D array
            self.flashes = np.row_stack((self.flashes, dataLine))
        print(self.flashes)
        return

    def _calcCrossCorr(self, dataA, dataB, tA, tB, width=0.5):
        """
        This function calculates a cross correlation of the data with a
        data width in data points. 
        """
        corrIndA = np.where((dataA['dateTime'] <= tA + timedelta(seconds=width)) & 
                            (dataA['dateTime'] >= tA - timedelta(seconds=width)))[0]
        corrIndB = np.where((dataB['dateTime'] <= tB + timedelta(seconds=width)) & 
                            (dataB['dateTime'] >= tB - timedelta(seconds=width)))[0] 
        dosA = dataA['dos1rate'][corrIndA]
        dosB = dataB['dos1rate'][corrIndB]
        # Replace any error values with 0 so correlate does not add any weight to it.
        dosA[dosA == -1E31] = 0
        dosB[dosB == -1E31] = 0
        cc = scipy.signal.correlate(dosA - np.mean(dosA), dosB - np.mean(dosB), 
            mode='valid')[0]
        cc /= np.std(dosA)*np.std(dosB)*(len(corrIndA)+len(corrIndB))/2
        return cc

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
            self.keys[2:], ['cross_correlation']))

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
    version = 2
    startTime = time.time()
    inDir = os.path.abspath('./../data/microburst_catalogues/')
    outDir = os.path.abspath('./../data/{}_catalogues/'.format(burstType))

    if 'flash' in burstType:
        flashesFlag = True
    else:
        flashesFlag = False
    sorter = SortMicrobursts(inDir, outDir, flashesFlag=flashesFlag)
    sorter.simpleFindMatches()
    sorter.saveData(outName='{}_catalogue_v{}.txt'.format(burstType, version))
    print('Run time: {}'.format(round(time.time()-startTime)))
