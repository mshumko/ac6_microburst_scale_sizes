# This script implements John's request to make a scatter plot of AE vs AC-6 
# separation.

import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser
import os

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

    def saveData(self, save_path):
        """ Saves data using pandas """
        df = pd.concat(
                [self.dataDict, pd.DataFrame(self.matchingIndex, columns=['AE'])], 
                axis=1, names=self.iType)
        df.to_csv(save_path, index=False)
        return 

if __name__ == '__main__':
    iType = 'ae'
    indexDir = '/home/mike/research/geomag_indicies/ae'
    separation_load_file = '/home/mike/research/ac6/AC6_Separation.csv'
    separation_save_file = '/home/mike/research/ac6/AC6_Separation_v2.csv'

    if not os.path.exists(separation_save_file):
        appendObj = AppendIndexToSeparation(iType, separation_load_file, indexDir)
        appendObj.appendIndex()
        appendObj.saveData(separation_save_file)
    else:
        df = pd.read_csv(separation_save_file)#, converters={0:dateutil.parser.parse})
        # plt.scatter(np.abs(df['In-Track Separation [km]']), df['AE'])
        # plt.xlabel('Separation [km]'); plt.ylabel('AE'); plt.title('AC6 AE vs Separation')
        # plt.show()

        # Make a pcolormesh plot of AE vs s.
        s_bins = np.arange(0, 200, 10)
        ae_bins = np.arange(0, 1000, 10)
        H, _, _ = np.histogram2d(np.abs(df['In-Track Separation [km]']), df['AE'], 
                                bins=(s_bins, ae_bins))
        im = plt.pcolormesh(s_bins[:-1], ae_bins[:-1], H.T)
        plt.colorbar(im, label='Number of intervals')
        plt.xlabel('Separation [km]'); plt.ylabel('AE'); plt.title('AC6 AE vs Separation')
        plt.show()

