# This code visualizes the MCMC trace models.

import numpy as np
import matplotlib.pyplot as plt
import os

import pandas as pd

import mcmc_models

TRACE_DIR = '/home/mike/research/ac6_microburst_scale_sizes/models/mcmc_traces'
TRACE_NAMES = ['analytic_trace_single_r.csv', 'mc_trace_single_r.csv']
CDF_DATA_PATH = ('/home/mike/research/ac6_microburst_scale_sizes'
            '/data/microburst_cdf_pdf_norm_v3.csv')

SIZE_BINS = np.arange(0, 100, 5)
cdf_data = pd.read_csv(CDF_DATA_PATH)

# Load the traces.
traces = {}
trace_stats = {}
# trace_stats = pd.DataFrame(data=np.nan*np.zeros((len(TRACE_NAMES), 3)), 
#                            columns=['mean', '2.5', '97.5'])
for i, trace_name in enumerate(TRACE_NAMES):
    traces[i] = pd.read_csv(os.path.join(TRACE_DIR, trace_name))

    # for key in traces[i]:
    #     trace_stats[i] = {}
    trace_stats[i] = traces[i].quantile([0.025, 0.5, 0.975])


# Visualize the traces.
fig, ax = plt.subplots(len(TRACE_NAMES), 2)
for i, trace_name in enumerate(TRACE_NAMES):
    ax[i, 0].hist(traces[i].to_numpy(), bins=SIZE_BINS)
    # ax[1, 1].plot()

plt.show()