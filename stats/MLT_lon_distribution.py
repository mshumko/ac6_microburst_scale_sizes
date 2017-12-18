# This script will plot the MLT-longitude distribution of flashes and curtains.

import os
import matplotlib.pyplot as plt

import global_distribution

Lrange=[4,6]

fig, axArr = plt.subplots(2, figsize=(10,8), sharex=True)
fPath = os.path.abspath('./../data/flash_catalogues/flashes_catalogue.txt')
pltObj = global_distribution.GlobalDistribution(fPath, burstType='flashes')
im = pltObj.rectLongDistribution(ax=axArr[0], Lrange=Lrange)
plt.colorbar(im, ax=axArr[0], label='Detections (unnormalized)')

fPath = os.path.abspath('./../data/curtain_catalogues/curtains_catalogue.txt')
pltObj = global_distribution.GlobalDistribution(fPath, burstType='curtains')
im2 = pltObj.rectLongDistribution(ax=axArr[1], Lrange=Lrange)
plt.colorbar(im2, ax=axArr[1], label='Detections (unnormalized)')

axArr[0].set(title='Flash distribution for {} < L < {}'.format(
    Lrange[0], Lrange[1]), ylabel='MLT')
axArr[1].set(title='Curtain distribution for {} < L < {}'.format(
    Lrange[0], Lrange[1]), ylabel='MLT', xlabel=r'Longitude [Deg]')

plt.tight_layout()
plt.show()