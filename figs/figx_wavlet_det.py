# This script generates a figure that shows how the wavelet-based microburst
# detector code works.

import numpy as np
import sys
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.ticker

sys.path.insert(0, '/home/mike/research/ac6-microburst-scale-sizes'
                    '/microburst-detection')
import microburst_detection

### DETECT MICROBURSTS ###
thresh=0.1
date = datetime(2016, 9, 30) 
obj = microburst_detection.FindMicrobursts('a', date)
obj.getMicroburstIdx(thresh=thresh)

### PLOT DETECTiONS ###
fig, ax = plt.subplots(3, sharex=True)
tRange = (datetime(2016, 9, 30, 0, 4, 51), datetime(2016, 9, 30, 0, 5, 41))
validIdt = np.where((obj.d['dos1rate'] != -1E31) & 
                    (obj.d['dateTime'] > tRange[0]) & 
                    (obj.d['dateTime'] < tRange[1]))[0]
validData = np.where(obj.d['dos1rate'] != -1E31)[0]
# Panel (a) shows the dos1rate points where detections were made.                   
ax[0].plot(obj.d['dateTime'][validIdt], 
                obj.d['dos1rate'][validIdt], label='dos1')               
#ax[0].fill_between(obj.d['dateTime'][validIdt], 
#    obj.d['dos1rate'][validIdt]-np.sqrt(obj.d['dos1rate'][validIdt]),
#    obj.d['dos1rate'][validIdt]+np.sqrt(obj.d['dos1rate'][validIdt]),
#   color='r', alpha=0.5)
ax[0].scatter(obj.d['dateTime'][validData[obj.burstIdt]], obj.d['dos1rate'][validData[obj.burstIdt]], c='b', s=50)
ax[0].scatter(obj.d['dateTime'][obj.peakInd], obj.d['dos1rate'][obj.peakInd],
    c='r', s=25)
#obj.plotPower(ax=ax[1])

### PLOT WAVELET ###
pIdx = np.where(obj.period < 3)[0]
levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
CS = ax[1].contourf(obj.time, obj.period[pIdx], np.log2(obj.power[pIdx, :]), len(levels))
ax[1].fill_between(obj.time, np.max(obj.period[pIdx]), obj.coi, alpha=0.8, 
    facecolor='white', zorder=3)
im = ax[1].contourf(CS, levels=np.log2(levels))
ax[1].contour(obj.time, obj.period[pIdx], obj.sig95[pIdx, :], [-99, 1], colors='k')
ax[1].axhspan(1, 14, alpha=0.7, color='grey')


ax[1].set_yscale('log', basey=2, subsy=None)
ax[1].set_ylim([np.min(obj.period[pIdx]), np.max(obj.period[pIdx])])
axy = ax[1].yaxis # plt.gca().yaxis
axy.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1].ticklabel_format(axis='y', style='plain')
ax[1].invert_yaxis()



ax[2].plot(obj.time, obj.dataFlt, label='filtered time series')
ax[2].axhline(thresh, c='k', label='threshold')


#self.ax[0].scatter(self.d['dateTime'][validIdt[self.burstIdt]], self.d['dos1rate'][validIdt[self.burstIdt]], c='b', s=50)
#self.ax[0].scatter(self.d['dateTime'][self.peakInd], self.d['dos1rate'][self.peakInd], c='r', s=25)


#obj.saveData()
#obj.plotTimeseries()
ax[0].set_xlim(*tRange)
ax[0].set_ylim(0.9*np.min(obj.d['dos1rate'][validIdt]), 
        1.1*np.max(obj.d['dos1rate'][validIdt]))
ax[2].set_ylim(-0.25, 0.5)        

fig.autofmt_xdate()
plt.tight_layout()
plt.show()                    
