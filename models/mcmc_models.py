# This program contains a few different microburst Monte Carlo (MC) models.

import numpy as np
import matplotlib.pyplot as plt
#import scipy.optimize
import scipy.stats
#import scipy.special
#import time
import os

# import theano.tensor as tt
# import pymc3 as pm
import pandas as pd 
import progressbar

GRID_SIZE = 200
# csv save data path. Will NOT overwrite if it already exists!
SAVE_PATH = ('/home/mike/research/ac6_microburst_scale_sizes/models/mcmc_traces'
            '/mc_trace_two_r.csv')
CDF_DATA_PATH = ('/home/mike/research/ac6_microburst_scale_sizes'
            '/data/microburst_cdf_pdf_norm_v3.csv')

# Define the analytic models.
def A(r, d): 
    """ 
    Calculates the intersecting area between two cirlces from the 
    analytic_coincident_microburst_scale_size_bounds.svg file.

    Derivation: http://mathworld.wolfram.com/Circle-CircleIntersection.html
    """
    return 2*r**2*np.arccos(d/(2*r)) - d/2*np.sqrt(4*r**2 - d**2)

def F(r, d):
    """ 
    Wrapper for A(r, d). Thus function calculates the fraction of coincident 
    microbursts observed at separation array values d. The cdf valie is 1 at 
    d[0] and monotonically approaches 0.
    """
    cdf = np.array([np.nansum(A(r, d[i:]))/np.nansum(A(r, d)) 
                    for i in range(len(d))])
    return cdf 

def distance(x1, y1, x2, y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)

def mc_brute_force(r, n_bursts=int(1E5), bins=np.arange(0, 100, 5)):
    """ 
    Brute force MC method that is computationally slow
    but should be the closest to reality.
    """
    N = np.zeros_like(bins)
    
    # If supplied a single valued radius r.
    if not hasattr(r, '__len__'):
        r = r*np.ones(n_bursts)
    
    for i, bin_i in enumerate(bins):
        # Generate n_bursts number of microbursts randomly scattered in a grid
        # and fixed radius r.
        burst_x = np.random.uniform(-GRID_SIZE, GRID_SIZE, size=n_bursts)
        burst_y = np.random.uniform(-GRID_SIZE, GRID_SIZE, size=n_bursts)
        
        # Now loop over the bursts and tally up the number of microbursts 
        # observed by hypothetical spacercaft at (0, 0) and (0, bin_i).
        for bx, by, br in zip(burst_x, burst_y, r):
            if (distance(bx, by, 0, 0) <= br) and (distance(bx, by, 0, bin_i) <= br): 
                N[i] += 1
    total_N = np.sum(N)
    cdf = np.array([np.sum(N[i:])/total_N for i in range(len(bins))])
    return cdf

def metroplis(start, target, proposal, niter, nburn=0, 
            thin=1, verbose=False):
    """
    This function implements the Metropolis–Hastings sampler.
    """
    niter = int(niter)
    current = start
    post = -1E31*np.ones((niter, len(start)), dtype=float)

    for i in progressbar.progressbar(range(niter)):
        proposed = proposal(current)
        p = min(target(proposed)/target(current), 1)
        if np.random.random() < p:
            current = proposed
        post[i, :] = current
        if verbose:
            print(f'Current values {current}')
    return post[int(nburn)::thin, :]

def Likelihood(p, x, y):
    """ Gaussian likelihood. """
    C = (np.std(y)*np.sqrt(2*np.pi))
    y_model = (p[0]*mc_brute_force(p[1], bins=x) + 
              (1-p[0])*mc_brute_force(p[2], bins=x))
    args = sum([(y_i - y_model_i)**2 
                for (y_i, y_model_i) in zip(y, y_model)])
    return np.exp(-0.5*args/np.var(y))/C


if __name__ == '__main__':
    if os.path.exists(SAVE_PATH):
        raise ValueError('Data already saved. Use a different savePath. Aborting.')

    # Load the CDF data to model
    cdf_data = pd.read_csv(CDF_DATA_PATH)
        
    # Specify priors for the one fixed-sized microburst model.
    # Two parameter model. First parameter is the mixing term, and second and third 
    # parameters are the two microburst sizes.
    prior = [scipy.stats.uniform(0, 1), 
            scipy.stats.uniform(0, 100), 
            scipy.stats.uniform(0, 100)]
    # Initial guess on the microburst size.
    start = [prior_i.rvs() for prior_i in prior]
    # How much to jump. Assuming a N(mu, sigma) proposal, ~60% of the time the next proposed jump will be less than proposal_jump km away.
    proposal_jump = [0.1, 2, 2]

    def proposal(p, proposal_jump=proposal_jump):
        new_vals = [scipy.stats.norm(loc=p_i, scale=jump_i).rvs() 
                        for p_i, jump_i in zip(p, proposal_jump)]
        # If the mixing term is not between 0 and 1, force it to 1 or 0.                  
        if new_vals[0] > 1:
            new_vals[0] = 1
        elif new_vals[0] < 0:
            new_vals[0] = 0
        return new_vals

    # The target function. If probability is higher, take the new value given from proposal. Else do the Metroplis thing where you draw a random number between 
    # 0 and 1 and compare to the target value (which will be less than 1).
    target = lambda p: Likelihood(p, cdf_data['Separation [km]'], 
                                cdf_data['CDF'])*np.prod(
                                [prior_i.pdf(p_i) for prior_i, p_i in zip(prior, p)])
    niter = 100
    trace = metroplis(start, target, proposal, niter, 
                            nburn=10, thin=1, verbose=False)
    # Save data
    df = pd.DataFrame(data=trace, columns=['a', 'r0', 'r1'])
    df.to_csv(SAVE_PATH, index=False)