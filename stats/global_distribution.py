# This script will plot the L-MLT distribution of flashes and curtains
import os
import numpy as np
import matplotlib.pyplot as plt

import scale_sizes

class GlobalDistribution(scale_sizes.ScaleSize):
    def __init__(self, fPath, burstType='flashes'):
        scale_sizes.ScaleSize.__init__(self, fPath, burstType)
        # self.burstType is set in the ScaleSize class.
        return

    def rectDistributionPlot(self, ax=None, Lbins=range(2, 10), MLTbins=range(24)):
        """ 
        This function plots the L-MLT distribution of flashes
        and curtains in a rectanguar plot.
        """
        if ax is None:
            _, self.ax = plt.subplots()
        else:
            self.ax = ax

        (counts, xedges, yedges, Image) = self.ax.hist2d(
            self.dataDict['MLT_OPQ'], self.dataDict['Lm_OPQ'], bins=[MLTbins, Lbins])
        if ax is None:
            self.ax.set(xlabel='MLT', ylabel='L', 
                title='L-MLT distribution of {}'.format(self.burstType))
            plt.colorbar(Image)
            plt.show()
        return Image

    def polarDistributionPlot(self, ax=None, Lbins=range(2, 10), MLTbins=range(24)):
        """
        This function plots the L-MLT distribution of flashes
        and curtains in a polar plot.
        """

        return 

if __name__ == '__main__':
    burstType = 'curtain'
    if 'flash' in burstType:
        fPath = os.path.abspath('./../data/flash_catalogues/flashes_catalogue.txt')
    else:
        fPath = os.path.abspath('./../data/curtain_catalogues/curtains_catalogue.txt')

    pltObj = GlobalDistribution(fPath, burstType=burstType+'s')
    #pltObj.rectDistributionPlot()