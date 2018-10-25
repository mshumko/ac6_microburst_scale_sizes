# These functions test the sort_curtains_flashes.py program
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np
import os
import sys

sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import sort_curtains_flashes
import read_ac_data

inDir = os.path.abspath('./../data/microburst_catalogues/')
date = datetime(2015, 11, 12)
sorter = sort_curtains_flashes.SortMicrobursts(inDir, date)

# Load AC6 data
dA = read_ac_data.read_ac_data_wrapper('A', date,
            dType='10Hz', plot=False)
dB = read_ac_data.read_ac_data_wrapper('B', date,
            dType='10Hz', plot=False)

# Plot the times to understand how closely they match
fig, ax = plt.subplots(4, sharex=True)

# Plot AC6 data
# Flashes
validIda = np.where(dA['dos1rate'] != -1E31)[0]
validIdb = np.where(dB['dos1rate'] != -1E31)[0]
ax[0].plot(dA['dateTime'][validIda], dA['dos1rate'][validIda], c='r', label='AC6A')
ax[0].plot(dB['dateTime'][validIdb], dB['dos1rate'][validIdb], c='b', label='AC6B')
ax[0].legend()

# Curtains
validIda = np.where(dA['dos1rate'] != -1E31)[0]
validIdb = np.where(dB['dos1rate'] != -1E31)[0]
# Shift times
dA['dateTime_shifted'] = np.array([dA['dateTime'][i] + timedelta(
    seconds=dA['Lag_In_Track'][i]) for i in validIda])
ax[2].plot(dA['dateTime_shifted'], dA['dos1rate'][validIda], c='r', label='AC6A')
ax[2].plot(dB['dateTime'][validIdb], dB['dos1rate'][validIdb], c='b', label='AC6B')

for burst in sorter.dataA['dateTime']:
    ax[1].axvline(burst, ymax=0.5, c='r')
for burst in sorter.dataB['dateTime']:
    ax[1].axvline(burst, ymin=0.5, ymax=1, c='b')

for i, burst in enumerate(sorter.dataA['dateTime']):
    ax[3].axvline(burst+timedelta(seconds=sorter.dataA['Lag_In_Track'][i]),
        ymax=0.5, c='r')
for burst in sorter.dataB['dateTime']:
    ax[3].axvline(burst, ymin=0.5, ymax=1, c='b')

ax[0].set_xlim(sorter.dataA['dateTime'][0], sorter.dataA['dateTime'][-1])
ax[0].set_title('Flashes')
ax[1].set_title('Curtains')
plt.show()
