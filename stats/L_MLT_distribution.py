# This code wrapps global_distribution.py and plots the L-MLT distributution
# for curtains and flashes.

import os
import matplotlib.pyplot as plt

import global_distribution

fig, axArr = plt.subplots(2, figsize=(15,8), sharex=True)
fPath = os.path.abspath('./../data/flash_catalogues/flashes_catalogue.txt')
pltObj = global_distribution.GlobalDistribution(fPath, burstType='flashes')
im = pltObj.rectDistributionPlot(ax=axArr[0])
plt.colorbar(im, ax=axArr[0], label='Detections (un-normalized)')

fPath = os.path.abspath('./../data/curtain_catalogues/curtains_catalogue.txt')
pltObj = global_distribution.GlobalDistribution(fPath, burstType='curtains')
im2 = pltObj.rectDistributionPlot(ax=axArr[1])
plt.colorbar(im2, ax=axArr[1], label='Detections (un-normalized)')

axArr[0].set(title='Flash distribution', ylabel='L')
axArr[1].set(title='Curtain distribution', ylabel='L', xlabel='MLT')

for ax in axArr:
    ax.set_aspect('equal')
plt.tight_layout()
plt.show()