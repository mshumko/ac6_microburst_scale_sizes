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

    def rectDistributionPlot(self, ax=None, Lbins=range(2, 10), MLTbins=range(24),
                            cmin=None, cmax=None):
        """ 
        This function plots the L-MLT distribution of flashes
        and curtains in a rectanguar plot.
        """
        if ax is None:
            _, self.ax = plt.subplots()
        else:
            self.ax = ax

        (counts, xedges, yedges, Image) = self.ax.hist2d(
            self.dataDict['MLT_OPQ'], self.dataDict['Lm_OPQ'], 
            bins=[MLTbins, Lbins], cmin=cmin, cmax=cmax)
        if ax is None:
            self.ax.set(xlabel='MLT', ylabel='L', 
                title='L-MLT distribution of {}'.format(self.burstType))
            plt.colorbar(Image)
            plt.show()
        return Image

    def polarDistributionPlot(self, ax=None, Lbins=range(3, 10), MLTbins=range(24)):
        """
        This function plots the L-MLT distribution of flashes
        and curtains in a polar plot.
        """
        (H, MLTedges, Ledges) = np.histogram2d(self.dataDict['MLT_OPQ'], self.dataDict['Lm_OPQ'],
            bins=[MLTbins, Lbins])

        width = (2*np.pi)/len(MLTbins)
        if ax is None:
            self.ax = plt.subplot(111, polar=True)
        else:
            self.ax = ax

        print(MLTedges, Ledges)
        bars = self.ax.bar(MLTedges, height=24, width=width, bottom=Ledges[0:-1])
        # Use custom colors and opacity
        for (ii, bar) in zip(H, bars):
            bar.set_facecolor(plt.cm.plasma(ii))
            #bar.set_alpha(0.8)

        plt.show()
        return 

if __name__ == '__main__':
    # burstType = 'flashes'
    # if 'flash' in burstType:
    #     fPath = os.path.abspath('./../data/flash_catalogues/flashes_catalogue.txt')
    # else:
    #     fPath = os.path.abspath('./../data/curtain_catalogues/curtains_catalogue.txt')

    #pltObj = GlobalDistribution(fPath, burstType=burstType)

    fig, axArr = plt.subplots(2, figsize=(15,8), sharex=True)
    fPath = os.path.abspath('./../data/flash_catalogues/flashes_catalogue.txt')
    pltObj = GlobalDistribution(fPath, burstType='flashes')
    im = pltObj.rectDistributionPlot(ax=axArr[0])
    plt.colorbar(im, ax=axArr[0], label='Detections (unformalized)')

    fPath = os.path.abspath('./../data/curtain_catalogues/curtains_catalogue.txt')
    pltObj = GlobalDistribution(fPath, burstType='curtains')
    im2 = pltObj.rectDistributionPlot(ax=axArr[1])
    plt.colorbar(im2, ax=axArr[1], label='Detections (unformalized)')

    axArr[0].set(title='Flash distribution', ylabel='L')
    axArr[1].set(title='Curtain distribution', ylabel='L', xlabel='MLT')


    #pltObj.rectDistributionPlot(ax=ax, cmin=0, cmax=120)
    # ax.set(xlabel='MLT', ylabel='L', 
    #             title='L-MLT distribution of {}'.format(burstType))
    for ax in axArr:
        ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()