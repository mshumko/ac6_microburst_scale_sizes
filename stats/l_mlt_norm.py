# This file will calculate the number of seconds AC6-B spent at each L-MLT
# bin

import sys
from datetime import datetime, timedelta
import itertools
import time
import csv
import numpy as np

import matplotlib.pylab as plt

sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

progStartTime = time.time()

startDate = datetime(2014, 1, 1)
endDate = datetime.now()
dDays = (endDate - startDate).days
dates = [startDate + timedelta(days=i) for i in range(dDays)]

sc_id = 'B'

# Create histogram arrays.
Lrange = np.arange(2, 11)
MLTrange = np.arange(0, 25)
timeHist = np.ones((len(Lrange)-1, len(MLTrange)-1), dtype=float)

# Loop over spcecraft and dates using itertools.product()
for date in dates: 
    print('Analyzing {}.'.format(date.date()))
    try:   
        data = read_ac_data.read_ac_data_wrapper(sc_id, date,
            dType='10Hz', plot=False)
    except AssertionError as err:
        if ( ('None or > 1 AC6 files found' in str(err)) or
            ('Error, the data is not 2D (Empty file?)' in str(err)) ):
            continue
        else:
            raise

    # Bin the data
    H, Ledges, MLTedges = np.histogram2d(data['Lm_OPQ'], data['MLT_OPQ'], 
        bins=[Lrange, MLTrange]) # Each data point is 0.1 s.
    timeHist += H/10
    
np.save('./normalization/L_MLT_sec_norm', timeHist)
np.save('./normalization/L_edges_norm', Ledges)
np.save('./normalization/MLT_edges_norm', MLTedges)

plt.imshow(timeHist)

plt.show()
