import scipy.stats
import scipy.integrate
import numpy as np

import pandas as pd

import mcmc_one_size


CDF_DATA_PATH = ('/home/mike/research/ac6_microburst_scale_sizes'
            '/data/microburst_cdf_pdf_norm_v3.csv')
cdf_data = pd.read_csv(CDF_DATA_PATH)

lower_bound=0
upper_bound=200
prior = [scipy.stats.uniform(lower_bound, upper_bound)]

def f(p):
    """ 
    Evaluate the likelihood*prior function and return a 
    non-zero value if not nan. 
    """
    L = mcmc_one_size.Likelihood([p],  
                                cdf_data['Separation [km]'], 
                                cdf_data['CDF'])
    x = L*prior[0].pdf(p)
    if np.isnan(x):
        return 0 
    else: 
        return x

def monte_carlo_integrator(f, priors, N):
    """ 
    This integrator plots over data points randomly chosen from the prior
    """
    # Generate N sets of prior paramers.
    p = np.nan*np.zeros((N, len(priors)))
    for i, prior in enumerate(priors):
        p[:, i] = prior.rvs(N)
    priors_bound = np.array([prior.interval(1) for prior in priors])
    V = np.prod(priors_bound[:, 1] - priors_bound[:, 0])
    f_eval = [f(p_i[0]) for p_i in p]
    integral = V/N*np.sum(f_eval)
    std = V*np.std(f_eval)/np.sqrt(N)
    return integral, std

integrand = [f(p_i) for p_i in np.arange(lower_bound, upper_bound+1)]
evidence = scipy.integrate.trapz(integrand, 
            x=np.arange(lower_bound, upper_bound+1))

mc_evidence, mc_std = monte_carlo_integrator(f, prior, int(1E3))

print(f'Trapz evidence = {evidence}')
print(f'MC evidence = {mc_evidence} +/- {mc_std}')