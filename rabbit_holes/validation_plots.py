# This script makes plots of the coincident micorbursts
import explore_ss_dependencies
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import os

plotPath = ('/home/mike/research/ac6-microburst-scale-sizes/'
            'data/plots/{}'.format(datetime.date.now()))
if not os.path.exists(plotPath):
    os.makedirs(plotPath)
    print('Made plot directory at', plotPath)




