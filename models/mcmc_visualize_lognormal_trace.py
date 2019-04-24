import numpy as np
import matplotlib.pyplot as plt
import os

import pandas as pd

TRACE_DIR = ('/home/mike/research/ac6_microburst_scale_sizes/models/'
             'mcmc_traces')
TRACE_NAME = 'mc_trace_log_norm.csv'
TRACE_PATH = os.path.join(TRACE_DIR, TRACE_NAME)

df = pd.read_csv(TRACE_PATH)

fig, ax = plt.subplots(2)
ax[0].plot(df.mu)
ax[1].plot(df.sigma)
plt.show()
