# This script finds the 10 most common bins for the microburst CDF as a function of separation.
import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser

import pandas as pd

# Load coincident microburst catalog
version = 4
catPath = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
        'coincident_microbursts_catalogues/'
       'AC6_coincident_microbursts_v{}.txt'.format(version))
converters = {
            0:dateutil.parser.parse, 
            -1:dateutil.parser.parse, 
            -2:dateutil.parser.parse
            }
data = pd.read_csv(catPath, converters=converters)

### FILTERS ###
# Filter out detections outside of the outer radiation belt
# and above the US and SAA. Also filter out detections
# which had a higher spatial_CC by a curtain_thresh ammount. 
# Lastly, filter out detections that are close to the noise
# floor (peak_std)
curtain_thresh = 0.3
# L filter
data = data[(data['Lm_OPQ'] > 4) & (data['Lm_OPQ'] < 8)]
# USA filter
data = data[
    ((data['lon'] > -60) | (data['lon'] < -140)) |
    ((data['lat'] > 70) | (data['lat'] < 15))
    ]
# SAA filter
data = data[
    ((data['lon'] > 30)  | (data['lon'] < -116)) |
    ((data['lat'] < -90) | (data['lat'] > 0))
    ]
# Filter by the number of standrad deviations a peak is 
# above a 10% baseline.
data = data[data['peak_std'] > 2]
# Filter out ambigious CCs with curtains
data = data[data['time_cc'] > data['space_cc']+curtain_thresh]
#bins = np.arange(2, 100, 3)