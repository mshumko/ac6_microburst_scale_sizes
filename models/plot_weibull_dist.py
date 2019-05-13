import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

PRIOR = 'norm'

# csv save data path. Will NOT overwrite if it already exists!
TRACE_PATH = ('/home/mike/research/ac6_microburst_scale_sizes/models/mcmc_traces'
            f'/mcmc_weibull_{PRIOR}_trace.csv')
df = pd.read_csv(TRACE_PATH)
#df = df[df.k < 2]

N = 5000
rand_ind = np.random.choice(df.shape[0], size=N)
s = np.linspace(0, 100, num=50)
pdf = np.nan*np.zeros((N, len(s)))

j = 0 # Counter

for _, row in df.iloc[rand_ind, :].iterrows():
    # parameter convention: c = k, scale=lambda, loc=offset
    w = scipy.stats.weibull_min(c=row.k, loc=row.offset, scale=row['lambda'])
    pdf[j, :] = w.pdf(s)
    #plt.plot(s, pdf[j, :]/max(pdf[j, :]), c='k', alpha=0.2)
    j += 1

percentiles = [50]
percentile_curves = np.percentile(pdf, percentiles, axis=0)
for i, percentile in enumerate(percentiles):
    plt.plot(s, percentile_curves[i, :])
#plt.legend()
plt.title(f'Median of {N} Weibull PDFs | MCMC with norm prior')
plt.xlabel('Microburst diamater [km]')
plt.ylabel('Microburst PD')
plt.show()