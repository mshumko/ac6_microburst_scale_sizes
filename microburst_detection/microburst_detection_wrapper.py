from datetime import datetime, timedelta, date
import itertools
import logging
import time
import os

import microburst_detection
import merge_daily_data
import replace_error_sep_lags
import append_goemag_indicies

cVersion = 5 #Catalog version
outPath = ('/home/mike/research/'
            'ac6-microburst-scale-sizes/data/'
            'microburst_catalogues/'
            'AC6A_microbursts_v{}.txt'.format(
            cVersion))
if os.path.exists(outPath):
    raise OSError('Merged catalog already exists. Increment cVersion.')

# Set up logger
logging.basicConfig(filename=os.path.abspath(
    './../logs/microburst_detection.log/'), level=logging.INFO, 
    format='%(asctime)s %(message)s')
    
progStartTime = time.time()

# Get a list of dates to loop over.
startDate = datetime(2014, 6, 21)
endDate = datetime(2017, 6, 30)
dDays = (endDate - startDate).days
dates = [startDate + timedelta(days=i) for i in range(dDays)]

# Loop over spcecraft and dates using itertools.product()
for (sc_id, date) in itertools.product(['A', 'B'], dates): 
   try:   
       # Set up detector and load AC-6 data.
       obj = microburst_detection.FindMicrobursts(sc_id, date)
       print('Analyzing {} on {}.'.format(sc_id, date.date()))
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
       # Run the wavelet detector.
#       obj.getMicroburstIdx(maxWidth=0.5, thresh=0.05,
#           SIGNIF_LEVEL=0.1) 
       obj.getMicroburstIdx(method='obrien')
       # Remove noisy detections using correlations.
       obj.corrFlag()
   except (ValueError, AssertionError) as err:
       if 'v cannot be empty' in str(err):
           logging.info('AC6-{} no microbursts found on {}.'.format(
               sc_id, date.date()))
           del(obj)
           continue
       elif 'attempt to get argmax of an empty sequence' in str(err):
           logging.info('AC6-{} no microbursts found on {}.'.format(
               sc_id, date.date()))
           del(obj)
           continue
       elif 'Data too short for timeWidth specified.' in str(err):
           logging.info('AC6-{} too little data for baseline. {}.'.format(
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
   inPath = ('/home/mike/research/ac6_microburst_scale_sizes/'
       'data/z_daily_microburst_catalogues')
       
   outPath = ('/home/mike/research/'
               'ac6_microburst_scale_sizes/data/'
               'microburst_catalogues/'
               'AC6{}_microbursts_v{}.txt'.format(
               sc_id.upper(), cVersion))
   if os.path.exists(outPath):
       raise OSError('Merged catalog already exists!')
   merge_daily_data.mergeDailyFiles(sc_id, inPath, outPath)

# Remove separation errors
for sc_id in ['A', 'B']:
    dName = 'AC6{}_microbursts_v{}.txt'.format(
            sc_id.upper(), cVersion)
    r = replace_error_sep_lags.ReplaceErrorVals(
        '/home/mike/research/ac6/AC6_Separation.csv', 
        ('/home/mike/research/ac6_microburst_scale_sizes/'
        'data/microburst_catalogues/{}'.format(dName)))
    r.loadSeprationFile()
    r.loadCatalog()
    r.replaceErrors()
    r.save_data()

# Append the AE index
iType = 'ae'
indexDir = '/home/mike/research/geomag_indicies/ae'
for sc_id in ['A', 'B']:
    dataPath = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
                'microburst_catalogues/'
                'AC6{}_microbursts_v{}.txt'.format(sc_id, cVersion))

    appendObj = append_goemag_indicies.AppendGeoMagIdx(iType, dataPath, indexDir)
    appendObj.appendIndex()
    appendObj.saveData()

logging.info('Microbust data pipeline ran in {}'.format(time.time() - progStartTime))
