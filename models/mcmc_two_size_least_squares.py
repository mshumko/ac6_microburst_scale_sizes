# Use the mcmc_one_size.py trace to estimate the model parameter 
# that minimize the least squares.

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import os

import pandas as pd 
import progressbar

import mcmc_two_size

# Load the trace and AC6 F(s)
PRIOR = 'uniform'
TRACE_PATH = ('/home/mike/research/ac6_microburst_scale_sizes/models/mcmc_traces'
            f'/mcmc_two_size_{PRIOR}_trace.csv')
CDF_DATA_PATH = ('/home/mike/research/ac6_microburst_scale_sizes'
            '/data/microburst_cdf_pdf_norm_v3.csv')
trace = pd.read_csv(TRACE_PATH)
cdf_data = pd.read_csv(CDF_DATA_PATH)

least_squares = np.nan*np.ones(trace.shape[0])
N_BURSTS = int(1E5)

#for i, p in progressbar.progressbar(trace.iterrows()):
for i in progressbar.progressbar(range(trace.shape[0])):
    p = trace.iloc[i]
    n_a = int(p.a*N_BURSTS)
    burst_diameters = np.concatenate((
        p.d0*np.ones(n_a),
        p.d1*np.ones(N_BURSTS-n_a)
        ))
    y_model = mcmc_two_size.mc_brute_vectorized(burst_diameters, 
                                bins=cdf_data['Separation [km]'], 
                                n_bursts=N_BURSTS)
    least_squares[i] = sum((cdf_data['CDF'] - y_model)**2)

optimal_trial = np.argmin(least_squares)
print(trace.iloc[optimal_trial])