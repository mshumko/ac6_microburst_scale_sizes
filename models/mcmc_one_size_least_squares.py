# Use the mcmc_one_size.py trace to estimate the model parameter 
# that minimize the least squares.

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import os

import pandas as pd 
import progressbar

import mcmc_one_size

# Load the trace and AC6 F(s)
PRIOR = 'uniform'
TRACE_PATH = ('/home/mike/research/ac6_microburst_scale_sizes/models/mcmc_traces'
            f'/mcmc_one_size_{PRIOR}_trace.csv')
CDF_DATA_PATH = ('/home/mike/research/ac6_microburst_scale_sizes'
            '/data/microburst_cdf_pdf_norm_v3.csv')
trace = pd.read_csv(TRACE_PATH)
cdf_data = pd.read_csv(CDF_DATA_PATH)

least_squares = np.nan*np.ones(trace.shape[0])

#for i, p in progressbar.progressbar(trace.iterrows()):
for i in progressbar.progressbar(range(trace.shape[0])):
    p = trace.iloc[i]
    y_model = mcmc_one_size.mc_brute_vectorized(p.d, bins=cdf_data['Separation [km]'])
    least_squares[i] = sum((cdf_data['CDF'] - y_model)**2)

optimal_trial = np.argmin(least_squares)
print(trace.iloc[optimal_trial])