"""
This script and functions looks at the microburst times in two
microburst catalogs and finds how many events are in both.
"""

import os
import numpy as np
import pandas as pd
from matplotlib.dates import date2num
from datetime import timedelta

cat_dir = ('/home/mike/research/ac6_microburst_scale_sizes/'
            'data/coincident_microbursts_catalogues')

def load_catalog(cat_name):
    """ 
    Load microburst catalog given by cat_name and 
    convert times to datetime objects 
    """
    df = pd.read_csv(os.path.join(cat_dir, cat_name))
    df['dateTime'] = pd.to_datetime(df['dateTime'])
    return df

def find_same_times(time_A, time_B, thresh=0.2):
    """ 
    This function finds the number of times both in time_A 
    and time_B with some small threshold. 
    """
    idx = 0
    for t_A in time_A:
        for t_B in time_B:
            if ((t_A - timedelta(seconds=thresh) < t_B) and 
                (t_A + timedelta(seconds=thresh) > t_B)):

                idx+=1
    return idx


cat_name_A = 'AC6_coincident_microbursts_sorted_v6.txt'
cat_name_B = 'AC6_coincident_microbursts_2scc_v9_sorted.txt'
cat_A = load_catalog(cat_name_A)
cat_B = load_catalog(cat_name_B)

idx = find_same_times(cat_A['dateTime'], cat_B['dateTime'])

print(f'Number of shared microbursts {idx}\n'
    f'Number of elements in A {cat_A.shape[0]}\n'
    f'Number of elements in B {cat_B.shape[0]}')