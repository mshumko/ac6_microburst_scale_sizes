# This program will sort microbursts by curtains and flashes and output the
# data to the ./../data/<flash or curtain>_catalogues.
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import dateutil.parser 
import csv
import os

class SortMicrobursts:
    def __init__(self, indir, date):
        self.date = date
        self._loadMicroburstCatalogues(indir)
        return

    def findFlashes(self):

        return

    def _loadMicroburstCatalogues(self, fDir):
        """
        This function will load in the microburst catalogues
        from both spacecraft
        """
        
        for sc_id in ['A', 'B']:
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
    inDir = os.path.abspath('./../data/microburst_catalogues/')
    date = datetime(2016, 9, 30)
    sorter = SortMicrobursts(inDir, date)

    # Plot the times to understand how closely they match
    fig, ax = plt.subplots()
    for burst in sorter.dataA['dateTime']:
        ax.axvline(burst, ymax=0.5, c='k')
    for burst in sorter.dataB['dateTime']:
        ax.axvline(burst, ymin=0.5, ymax=1, c='b')
    ax.set_xlim(sorter.dataA['dateTime'][0], sorter.dataA['dateTime'][-1])
    plt.show()