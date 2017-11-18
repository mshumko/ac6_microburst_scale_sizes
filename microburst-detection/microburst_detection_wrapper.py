from datetime import datetime, timedelta, date
import itertools
import logging
import time

import microburst_detection

# Set up logger
logging.basicConfig(filename=os.path.abspath(
    './../logs/microburst_detection.log/'), level=logging.INFO, 
    format='%(asctime)s %(message)s')
    
progStartTime = time.time()

startDate = date(2014, 1, 1)
endDate = datetime.now()
dDays = (endDate - startDate).days

dates = [startDate + timedelta(days=i) for i in range(dDays)]
# Loop over spcecraft and dates using itertools.product()
for (sc_id, date) in itertools.product(['A', 'B'], dates): 
    try:   
        obj = microburst_detection.TestFindMicrobursts(sc_id, date)
    except AssertionError as err:
        if 'None or > 1 AC6 files found' in err:
            logging.debug(err) # Wont log is logging level is logging.DEBUG
            continue
    obj.getMicroburstIdx()
    obj.saveData()
    logging.info('AC6-{}, {} microbursts detected'.format(sc_id, date))
    
logging.info('Program ran in {}'.format(time.time() - progStartTime))
