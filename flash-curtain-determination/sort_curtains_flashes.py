# This program will sort microbursts by curtains and flashes and output the
# data to the ./../data/<flash or curtain>_catalogues.
from datetime import datetime, timedelta
import numpy as np
import csv
import os

class SortMicrobursts:
    def __init__(self, indir, date):
        self.date = date
        self._loadMicroburstCatalogues(indir)
        return

    def _loadMicroburstCatalogues(self, fDir):
        """
        This function will load in the microburst catalogues
        from both spacecraft
        """
        
        for sc_id in ['A', 'B']:
            fName = 'AC6{}_{}_micobursts.txt'.format(
                sc_id, self.date.date())
            self._dataSkelletons(sc_id, os.path.join(fDir, fName))

        return

    def _dataSkelletons(self, sc_id, fPath):
        """
        Create empty dictionary and arrays for microburst catalogues. 
        """
        # Scan the file to get number of lines.
        with open(fPath, 'r') as f:
            N = sum(1 for row in f) - 2

        # Get data keys from file
        with open(fPath, 'r') as f:
            reader = csv.reader(f)
            next(reader) # Skip first line of header
            self.keys = next(reader)
            self.keys[0] = self.keys[0][2:] # Remove comment and empty space chars

            
            data = np.nan*np.ones((N, len(self.keys)))

            for i, line in enumerate(reader): # Read rest of file
                data[i] = line

            print(data)
        return

if __name__ == '__main__':
    inDir = os.path.abspath('./../data/microburst_catalogues/')
    date = datetime(2016, 9, 30)
    sorter = SortMicrobursts(inDir, date)