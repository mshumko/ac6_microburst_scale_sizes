# This function will claculate the normalization coefficients for the 
# scale_sizes by the number of days in each separation bin.
import sys
from datetime import datetime, timedelta
import itertools
import time
import csv
import numpy as np

sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

progStartTime = time.time()

startDate = datetime(2014, 1, 1)
endDate = datetime.now()
dDays = (endDate - startDate).days

sc_id = 'B'
sep_bins = np.arange(0, 500, 5)
sep_count = np.zeros_like(sep_bins)

dates = [startDate + timedelta(days=i) for i in range(dDays)]
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
    
    # Now calculate the mean total separation, and bin it.
    validD = np.where(data['Dist_Total'] != -1E31)[0]
    dist = np.mean(data['Dist_Total'][validD])
    
    for i, (bL, bH) in enumerate(zip(sep_bins[:-1], sep_bins[1:])):
        if (dist > bL) and (dist < bH):
            sep_count[i] += 1
            break

with open('dist_norm.txt', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['# AC6 daily average separation'])
    writer.writerow(['# Dist_Total (km)', 'Count (days)'])
    
    for line in zip(sep_bins, sep_count):
        writer.writerow(line)
