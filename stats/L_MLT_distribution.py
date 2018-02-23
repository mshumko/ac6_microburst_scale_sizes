# This code wrapps global_distribution.py and plots the L-MLT distributution
# for curtains and flashes.

import os
import matplotlib.pyplot as plt
import numpy as np

import global_distribution

# Load the normalization histogram
L_edges = np.load(os.path.abspath('./normalization/L_edges_norm.npy'))
MLT_edges = np.load(os.path.abspath('./normalization/MLT_edges_norm.npy'))
norm = np.load(os.path.abspath('./normalization/L_MLT_sec_norm.npy')).T+1

fig, axArr = plt.subplots(2, figsize=(15,8), sharex=True)
fPath = os.path.abspath('./../data/flash_catalogues/flash_catalogue_v2_sorted.txt')
pltObj = global_distribution.GlobalDistribution(fPath, burstType='flashes')
im = pltObj.rectDistributionPlot(ax=axArr[0], Lbins=L_edges, MLTbins=MLT_edges,
    norm=norm)
plt.colorbar(im, ax=axArr[0], label='Detections (un-normalized)')

# fPath = os.path.abspath('./../data/curtain_catalogues/curtains_catalogue.txt')
# pltObj = global_distribution.GlobalDistribution(fPath, burstType='curtains')
# im2 = pltObj.rectDistributionPlot(ax=axArr[1], Lbins=L_edges, MLTbins=MLT_edges,
#     norm=norm)
# plt.colorbar(im2, ax=axArr[1], label='Detections (un-normalized)')

axArr[0].set(title='Flash distribution', ylabel='L')
#axArr[1].set(title='Curtain distribution', ylabel='L', xlabel='MLT')

for ax in axArr:
    ax.set_aspect('equal')
plt.tight_layout()
plt.show()