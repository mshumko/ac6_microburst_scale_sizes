import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import os

import pandas as pd 
import progressbar

import mcmc_models

# csv save data path. Will NOT overwrite if it already exists!
SAVE_PATH = ('/home/mike/research/ac6_microburst_scale_sizes/models/mcmc_traces'
            '/mc_trace_log_norm.csv')
CDF_DATA_PATH = ('/home/mike/research/ac6_microburst_scale_sizes'
            '/data/microburst_cdf_pdf_norm_v3.csv')

if os.path.exists(SAVE_PATH):
        raise ValueError('Data already saved. Use a different savePath. Aborting.')

# Load the CDF data to model
cdf_data = pd.read_csv(CDF_DATA_PATH)
    
# Specify priors for the log-normal distribution. First parameter is the mean and second
#  is the standard deviation
prior = [scipy.stats.uniform(0, 100), 
        scipy.stats.uniform(0, 10)]
# Initial guess on the microburst size.
start = [prior_i.rvs() for prior_i in prior]
# How much to jump. Assuming a N(mu, sigma) proposal, ~60% of the time the 
# next proposed jump will be less than proposal_jump km away.
proposal_jump = [2, 1]

def proposal(p, proposal_jump=proposal_jump):
    new_vals = [scipy.stats.norm(loc=p_i, scale=jump_i).rvs() 
                    for p_i, jump_i in zip(p, proposal_jump)]
    # If the mixing term is not between 0 and 1, force it to 1 or 0.                  
    return new_vals

def Likelihood(p, x, y):
    """ Gaussian likelihood. """
    r = np.random.lognormal(mean=np.log(p[0]), sigma=p[1], size=int(1E5))
    C = (np.std(y)*np.sqrt(2*np.pi))
    y_model = mcmc_models.mc_brute_force(2*r, bins=x)
    args = sum([(y_i - y_model_i)**2 
                for (y_i, y_model_i) in zip(y, y_model)])
    return np.exp(-0.5*args/np.var(y))/C

# The target function. If probability is higher, take the new value given from proposal. Else do the Metroplis thing where you draw a random number between 
# 0 and 1 and compare to the target value (which will be less than 1).
target = lambda p: Likelihood(p, cdf_data['Separation [km]'], 
                            cdf_data['CDF'])*np.prod(
                            [prior_i.pdf(p_i) for prior_i, p_i in zip(prior, p)])
niter = 1000
trace = mcmc_models.metroplis(start, target, proposal, niter, 
                        nburn=100, thin=1, verbose=False)
# Save data
df = pd.DataFrame(data=trace, columns=['mu', 'sigma'])
df.to_csv(SAVE_PATH, index=False)