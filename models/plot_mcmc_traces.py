# This code visualizes the MCMC trace models.

import numpy as np
import matplotlib.pyplot as plt
import os

import pandas as pd

import mcmc_models

TRACE_DIR = '/home/mike/research/ac6_microburst_scale_sizes/models/mcmc_traces'
TRACE_NAMES = ['analytic_trace_single_r.csv', 'mc_trace_single_r.csv'] #'mc_trace_two_r.csv']
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
fig, ax = plt.subplots(len(TRACE_NAMES), 2, figsize=(8, 8))
for i, trace_name in enumerate(TRACE_NAMES):
    ax[i, 0].hist(traces[i].to_numpy(), bins=SIZE_BINS)
    ax[i, 0].set_ylabel(trace_name.split('.')[0])

    # Plot the CDF data.
    ax[i, 1].plot(cdf_data['Separation [km]'], cdf_data['CDF'], c='k', label='_nolegend_')
    ax[i, 1].fill_between(cdf_data['Separation [km]'], 
                 cdf_data['CDF']-cdf_data['CDF_std'], 
                 cdf_data['CDF']+cdf_data['CDF_std'], 
                 color='k', alpha=0.3, label='data')

# Plot the fits to the data from the mean and 95% values.
for i in trace_stats[0].to_numpy():
    ax[0, 1].plot(SIZE_BINS, mcmc_models.F(i[0], SIZE_BINS))

for i in trace_stats[1].to_numpy():
    ax[1, 1].plot(SIZE_BINS, mcmc_models.mc_brute_force(i[0], bins=SIZE_BINS))

#for row in trace_stats[2].iterrows():
#    cdf = (row[1].a*mcmc_models.mc_brute_force(row[1].r0, bins=SIZE_BINS) + 
#          (1-row[1].a)*mcmc_models.mc_brute_force(row[1].r1, bins=SIZE_BINS))
#    ax[2, 1].plot(SIZE_BINS, cdf)
#    # row[1].r0

plt.tight_layout()
plt.show()
