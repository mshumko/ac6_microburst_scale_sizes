# This script makes plots of the coincident micorbursts
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import os
import sys


sys.path.insert(0, '/home/mike/research/mission-tools/ac6')

import read_ac_data

# Path containing the catalog file to validate
catPath = os.path.join('/home/mike/research/ac6-microburst-scale-sizes/data/'
        'coincident_microbursts_catalogues', 'flash_catalogue_v2_sorted.txt')
pltWidth = 5 # seconds

plotPath = ('/home/mike/research/ac6-microburst-scale-sizes/'
            'data/plots/{}'.format(datetime.date.now()))
if not os.path.exists(plotPath):
    os.makedirs(plotPath)
    print('Made plot directory at', plotPath)

# Load and filter the catalog    
c = ExploreDependencies(catPath)
c.filter()

# Set up plots
fig, ax = plt.subplots(1)
current_date = datetime.date().min
for i in len(c.cat['burstType']):
    if c.cat['dateTimeA'][i].date().isoformat != current_date:
        dataA = read_ac_data_wrapper('A', current_date, dType='10Hz')
        dataB = read_ac_data_wrapper('B', current_date, dType='10Hz')
        current_date = c.cat['dateTimeA'][i].date().isoformat
        
        # Pick out only the valid data
        validIdxA = np.where(dataA['dos1rate'] > 0)[0]
        validIdxB = np.where(dataB['dos1rate'] > 0)[0]
        # Plot the unshifted data
        ax[0].plot(dataA['dateTime'][validIdxA], dataA['dos1rate'][validIdxA], 
                    label='AC-6 A')
        ax[0].plot(dataB['dateTime'][validIdxB], dataB['dos1rate'][validIdxB], 
                    label='AC-6 B')    
                    
    # Set the time range around the coincident microburst event, and then save.
    xlim=(c.cat['dateTimeA'][i] - timedelta(seconds=pltWidth), 
            c.cat['dateTimeA'][i] + timedelta(seconds=pltWidth))
           
    # Calculate ylimits from the xlimits
    idxA = np.where((dataA['dateTime'] > xlim[0]) & 
                    (dataA['dateTime'] < xlim[1]) &
                    (dataA['dos1rate'] > 0))[0]
    idxB = np.where((dataB['dateTime'] > xlim[0]) & 
                    (dataB['dateTime'] < xlim[1]) &
                    (dataB['dos1rate'] > 0))[0]
    ylim = (0.9*np.min(dataA['dos1rate'][idxA]), 
            1.1*np.max(dataA['dos1rate'][idxA]))
            
    ax[0].set(title='AC-6 Coincident microburst validation', xlabel='UTC', 
            ylabel='dos1 [counts/s]', xlim=xlim, ylim=ylim)
    ax[0].legend(loc=1)
    
    plt.savefig(os.path.join(plotPath, ))
        
        
    
