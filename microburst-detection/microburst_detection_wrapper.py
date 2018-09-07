from datetime import datetime, timedelta, date
import itertools
import logging
import time
import os

import microburst_detection
import merge_daily_data

cVersion = 3 #Catalog version

# Set up logger
logging.basicConfig(filename=os.path.abspath(
    './../logs/microburst_detection.log/'), level=logging.INFO, 
    format='%(asctime)s %(message)s')
    
progStartTime = time.time()

# Get a list of dates to loop over.
startDate = datetime(2014, 1, 1)
endDate = datetime.now()
dDays = (endDate - startDate).days
dates = [startDate + timedelta(days=i) for i in range(dDays)]

# Loop over spcecraft and dates using itertools.product()
for (sc_id, date) in itertools.product(['A', 'B'], dates): 
    print('Analyzing {} on {}.'.format(sc_id, date.date()))
    try:   
        # Set up detector and load AC-6 data.
        obj = microburst_detection.FindMicrobursts(sc_id, date)
    except AssertionError as err:
        if ( ('None or > 1 AC6 files found' in str(err)) or
            ('Error, the data is not 2D (Empty file?)' in str(err)) 
            or ('File is empty!' in str(err))):
            # Wont log is logging level is logging.DEBUG
            logging.debug('AC6-{} on {}: {}'.format(sc_id, date.date(), err)) 
            continue
        else:
            raise
    try:
        obj.getMicroburstIdx(maxWidth=0.5) # Run the wavelet detector.
    except ValueError as err:
        if 'v cannot be empty' in str(err):
            logging.info('AC6-{} no microbursts found on {}.'.format(
                sc_id, date.date()))
            del(obj)
            continue
        if 'attempt to get argmax of an empty sequence' in str(err):
            logging.info('AC6-{} no microbursts found on {}.'.format(
                sc_id, date.date()))
            del(obj)
            continue
        else:
            raise
    obj.saveData()
    del(obj)
    logging.info('AC6-{}, {} microbursts detected'.format(sc_id, date.date()))
    
# Merge daily files
for sc_id in ['a', 'b']:
    inPath = ('/home/mike/research/ac6-microburst-scale-sizes/'
        'data/z_daily_microburst_catalogues')
        
    outPath = ('/home/mike/research/'
                'ac6-microburst-scale-sizes/data/'
                'microburst_catalogues/'
                'AC6{}_microbursts_v{}.txt'.format(
                sc_id.upper(), cVersion))
    if os.path.exists(outPath):
        raise OSError('Merged catalog already exists!')
    merge_daily_data.mergeDailyFiles(sc_id, inPath, outPath)

logging.info('Program ran in {}'.format(time.time() - progStartTime))
