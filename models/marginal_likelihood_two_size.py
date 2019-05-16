import scipy.stats
import scipy.integrate
import numpy as np

import pandas as pd

import mcmc_two_size


CDF_DATA_PATH = ('/home/mike/research/ac6_microburst_scale_sizes'
            '/data/microburst_cdf_pdf_norm_v3.csv')
cdf_data = pd.read_csv(CDF_DATA_PATH)

priors = [
        scipy.stats.uniform(0, 0.2), 
        scipy.stats.uniform(50, 200),
        scipy.stats.uniform(0, 50) 
        ]

# priors = [
#         scipy.stats.halfnorm(loc=0, scale=0.3), 
#         scipy.stats.norm(loc=100, scale=50),
#         scipy.stats.norm(loc=30, scale=20) 
#         ]

def f(p):
    """ 
    Evaluate the likelihood*prior function and return a 
    non-zero value if not nan. 
    """
    L = mcmc_two_size.gaus_likelihood(p,  
                                cdf_data['Separation [km]'], 
                                cdf_data['CDF'], int(1E5))
    x = L*np.prod([prior.pdf(p[i]) for i, prior in enumerate(priors)])
    if np.isnan(x):
        return 0 
    else: 
        # if np.inf == x:
        #     print(x)
        # print()
        #print(L, np.prod([prior.pdf(p[i]) for i, prior in enumerate(priors)]), x)
        return x

def monte_carlo_integrator(f, priors, N):
    """ 
    This integrator plots over data points randomly chosen from the prior
    """
    # Generate N sets of prior paramers.
    p = np.nan*np.zeros((N, len(priors)))
    for i, prior in enumerate(priors):
        p[:, i] = prior.rvs(N)
    # Fix the mixing term to 1 if it is greater than 1.
    idx = np.where(p[:, 0] > 1)[0]
    p[idx, 0] = 1
    # Calculate the integration volume.
    priors_bound = np.array([prior.interval(1) for prior in priors])
    V = np.prod(priors_bound[:, 1] - priors_bound[:, 0])

    # f_eval = np.nan*np.zeros(N)
    # for i in range(N):
    #     f_eval[i] = f(p[i, :])
    #     try:
            
    #         if f_eval[i] > 1:
    #             print(f_eval[i])
    #     except:
    #         print(p[i, :])
    #         raise
    f_eval = [f(p_i) for p_i in p]
    # for f_i in f_eval:
    #     print(f_i)
    # print(np.where(np.isinf(f_eval))[0])
    integral = V/N*np.sum(f_eval)
    std = V*np.std(f_eval)/np.sqrt(N)
    return f_eval, integral, std

# integrand = [f(p_i) for p_i in np.arange(lower_bound, upper_bound+1)]
# evidence = scipy.integrate.trapz(integrand, 
#             x=np.arange(lower_bound, upper_bound+1))

f_eval, mc_evidence, mc_std = monte_carlo_integrator(f, priors, int(1E3))

# print(f'Trapz evidence = {evidence}')
print(f'MC evidence = {mc_evidence} +/- {mc_std}')