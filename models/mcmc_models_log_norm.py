import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import os

import pandas as pd 
import progressbar

import mcmc_models

OVERWRITE = False

# csv save data path. Will NOT overwrite if it already exists!
SAVE_PATH = ('/home/mike/research/ac6_microburst_scale_sizes/models/mcmc_traces'
            '/mc_trace_log_norm.csv')
CDF_DATA_PATH = ('/home/mike/research/ac6_microburst_scale_sizes'
            '/data/microburst_cdf_pdf_norm_v3.csv')

# if os.path.exists(SAVE_PATH):
#         raise ValueError('Data already saved. Use a different savePath. Aborting.')

# Load the CDF data to model
cdf_data = pd.read_csv(CDF_DATA_PATH)
    
# Specify priors for the log-normal distribution. First parameter is the mean and second
#  is the standard deviation
prior = [scipy.stats.uniform(0, 100), 
        scipy.stats.uniform(0.1, 10)]
# Initial guess on the microburst size.
start = [prior_i.rvs() for prior_i in prior]
# How much to jump. Assuming a N(mu, sigma) proposal, ~60% of the time the 
# next proposed jump will be less than proposal_jump km away.
proposal_jump = [5, 0.1]

def proposal(p, proposal_jump=proposal_jump):
    new_vals = np.array([scipy.stats.norm(loc=p_i, scale=jump_i).rvs() 
                        for p_i, jump_i in zip(p, proposal_jump)])
    # If parameters are less than 0, force them to 0.    
    new_vals[new_vals < 0] = 0              
    return new_vals

def Likelihood(p, x, y):
    """ Gaussian likelihood. """
    r = np.random.lognormal(mean=np.log(p[0]), sigma=p[1], size=int(1E5))
    C = (np.std(y)*np.sqrt(2*np.pi))
    y_model = mcmc_models.mc_brute_vectorized(2*r, bins=x)
    args = sum([(y_i - y_model_i)**2 
                for (y_i, y_model_i) in zip(y, y_model)])
    return np.exp(-0.5*args/np.var(y))/C

# The target function. If probability is higher, take the new value given from proposal. 
# Else do the Metroplis thing where you draw a random number between 
# 0 and 1 and compare to the target value (which will be less than 1).
target = lambda p: Likelihood(p, cdf_data['Separation [km]'], 
                            cdf_data['CDF'])*np.prod(
                            [prior_i.pdf(p_i) for prior_i, p_i in zip(prior, p)])
niter = 100000

if OVERWRITE or (not os.path.exists(SAVE_PATH)):
    trace = mcmc_models.metroplis(start, target, proposal, niter, 
                            nburn=100, thin=1, verbose=False)
    # Save data
    df = pd.DataFrame(data=trace, columns=['mu', 'sigma'])
    df.to_csv(SAVE_PATH, index=False)
else:
    print('Data already saved. Aborting MCMC.')
    df = pd.read_csv(SAVE_PATH)

quantiles = [2.5, 50, 97.5]
print(df.quantile([0.025, 0.50, 0.975]))

### PLOTTING CODE ###
fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(3, 2)
ax = np.zeros((2, 2), dtype=object)
for row in range(ax.shape[0]):
    for column in range(ax.shape[1]):
        ax[row, column] = plt.subplot(gs[row, column])
bx = plt.subplot(gs[-1, :])

N = df.shape[0]
ax[0,0].plot(np.arange(N)/1E4, df.mu, c='k')
ax[1,0].hist(df.mu, density=True, bins=np.linspace(0, 100), color='k')
ax[0,0].set(xlabel=r'Iteration x $10^4$', ylabel='mu trace')
ax[1,0].plot(np.linspace(0, 100), prior[0].pdf(np.linspace(0, 100)))
ax[1,0].set(xlabel='mu', ylabel='mu posterior PD')

ax[0,1].plot(np.arange(N)/1E4, df.sigma, c='k')
ax[0,1].set(xlabel=r'Iteration x $10^4$', ylabel=r'sigma trace')
ax[1,1].hist(df.sigma, density=True, bins=np.linspace(0, 5), color='k')
ax[1,1].plot(np.linspace(0, 5), prior[1].pdf(np.linspace(0, 5)))
ax[1,1].set(xlabel=r'sigma', ylabel=r'sigma posterior PD')

# Pick 1000 traces to analyze further to make plots.
rand_ind = np.random.choice(np.arange(N), size=1000)
y_model = np.nan*np.zeros((len(rand_ind), 
                        len(cdf_data['Separation [km]'])))

N_plot = 100
j = 0
for _, row in df.iloc[rand_ind, :].iterrows():
    #print(row)
    burst_diameters = np.random.lognormal(mean=np.log(row.mu), 
                                        sigma=row.sigma, size=100000)
    y_model[j, :] = mcmc_models.mc_brute_vectorized(burst_diameters, 
                            bins=cdf_data['Separation [km]'])
    j += 1

for i in range(N_plot):
    bx.plot(cdf_data['Separation [km]'], y_model[i,:], c='grey', alpha=0.2)
    bx.plot(cdf_data['Separation [km]'], cdf_data['CDF'], c='k')

quantiles = [2.5, 50, 97.5]
y_quartile = np.nanpercentile(y_model, quantiles, axis=0)

colors = ['g', 'r', 'b']
for i, q in enumerate(y_quartile):
    bx.plot(cdf_data['Separation [km]'], q, c=colors[i], 
            label=f'{quantiles[i]}')

plt.suptitle('log-normal microburst population MCMC model')
bx.set(xlabel='Spacecraft separation [km]', ylabel='F(d)')
gs.tight_layout(fig)
plt.show()
