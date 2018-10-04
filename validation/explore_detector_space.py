# This script explores the parameter space of the wavelet and burst parameter
# detection algorithms to find which parameters highest true-positive rate 
# with lowest false-positive rates.  

import numpy as np
import csv
import confusion_matrix
from datetime import datetime
import itertools
import os

# Test wavelet params
thresh = np.arange(1, 20)/20
maxWidth = np.arange(0.3, 2, 0.1)
SIGNIF_LEVEL = 1/np.arange(1, 11) #[1/n for n in range(1, 10)]

c = confusion_matrix.ConfusionMatrix('A', datetime(2016, 10, 14))
save_name = 'wavelet_params.csv'
if not os.path.isfile(save_name):
    with open(save_name, 'w') as f:
        f.write('thresh, maxWidth, SIGNIF_LEVEL, validNum, detNum, TPR, FPR\n')


for (t, w, s) in itertools.product(thresh, maxWidth, SIGNIF_LEVEL):
    dKwargs = {'method':'wavelet', 'thresh':t, 'maxWidth':w,
                    'SIGNIF_LEVEL':s}
    #dKwargs = {'method':'obrien', 'n':0.1, 'a':0.5, 'thresh':5}
    c.find_microbursts(**dKwargs)
    c.confusionMatrix()
    
    with open(save_name, 'a') as f:
        f.write('{}, {}, {}, {}, {}, {}, {}\n'.format(
                round(t, 2), round(w, 1), round(s, 1), len(c.vNum), 
                len(c.peakInd), round(100*c.TPR), 
                round(100*c.numFalse/len(c.vNum))
                ))
    
