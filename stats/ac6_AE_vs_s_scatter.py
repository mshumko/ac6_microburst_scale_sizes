# This script implements John's request to make a scatter plot of AE vs AC-6 
# separation.

import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser

import pandas as pd

from ac6_microburst_scale_sizes.microburst_detection.append_goemag_indicies import AppendGeoMagIdx

class AppendIndexToSeparation(AppendGeoMagIdx):
    def __init__(self, iType, dataPath, indexDir):
        super().__init__(iType, dataPath, indexDir, timeKey='Date/Time')

    def _loadData(self):
        """ Overwrite method """
        self.dataDict = pd.read_csv(self.dataPath, 
                    converters={0:dateutil.parser.parse})
        self.keys = self.dataDict.keys() # For the parent class.
        return

    def saveData(self):
        """ Saves data using pandas """
        df = pd.concat(
                [self.dataDict, pd.DataFrame(self.matchingIndex, columns=['AE'])], 
                axis=1, names=self.iType)
        df.to_csv(self.dataPath, index=False)
        return 

if __name__ == '__main__':
    iType = 'ae'
    indexDir = '/home/mike/research/geomag_indicies/ae'
    separation_file = '/home/mike/research/ac6/AC6_Separation_v2.csv'

    appendObj = AppendIndexToSeparation(iType, separation_file, indexDir)
    appendObj.appendIndex()
    appendObj.saveData()


