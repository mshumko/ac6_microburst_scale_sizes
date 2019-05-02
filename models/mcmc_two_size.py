# This program run a MCMC simulation of a two fixed-sized microburst 
# populations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats
import os

import pandas as pd 
import progressbar

GRID_SIZE = 200
OVERWRITE = False
LIKELIHOOD_ERROR = 0.1

# csv save data path. Will NOT overwrite if it already exists!
SAVE_PATH = ('/home/mike/research/ac6_microburst_scale_sizes/models/mcmc_traces'
            '/mcmc_two_size_trace.csv')
CDF_DATA_PATH = ('/home/mike/research/ac6_microburst_scale_sizes'
            '/data/microburst_cdf_pdf_norm_v3.csv')
    
def mc_brute_vectorized(burst_diamaters, n_bursts=100000, 
                         bins=np.arange(0, 100, 5), 
                         grid_size=200):
    """ 
    This function tallies events that were see by both spacercaft over 
    the events seen by one, or both spacercaft as a function of spacecraft 
    separation.
    """
    n = np.zeros(len(bins))
    
    # If supplied a single-valued diameter burst_diamaters. If supplied an array, 
    # I implicityly assume the diameters are distributed according to 
    # some PDF distribution.
    if not hasattr(burst_diamaters, '__len__'):
        burst_diamaters = burst_diamaters*np.ones(n_bursts)
    # Randomly generate n_burst microburst centers.
    burst_centers = np.random.uniform(-grid_size, grid_size, size=(n_bursts, 2))
    
    ### Calculate which bursts interesect the origin. ###
    # First calculate the distance each burst is from the origin
    spacecraft_location = np.zeros_like(burst_centers)
    distance_to_origin = np.linalg.norm(burst_centers-spacecraft_location, axis=1)
    # Where close_to_origin is True, the microburst intersected the origin.
    # We only need to loop over the True values.
    close_to_origin = np.less(distance_to_origin, burst_diamaters/2)
    i_close_to_origin = np.where(close_to_origin)[0]
    # Filter the burst centers to loop over only the ones that contain the origin.
    burst_centers = burst_centers[i_close_to_origin, :]
    burst_diamaters = burst_diamaters[i_close_to_origin]
    
    # Loop over all spacecraft bins and calculate the subset of the microbursts
    # that contain the origin also contain the spacecraft at distance d from the
    # origin along one axis (positive y-axis in this model).
    for i, d in enumerate(bins):
        spacecraft_location = np.zeros_like(burst_centers)
        spacecraft_location[:, 1] = d
        distance_to_spacecraft = np.linalg.norm(burst_centers-spacecraft_location, 
                                            axis=1)
        # Where close_to_origin is True, the microburst intersected the origin.
        # We only need to loop over the True values.
        close_to_spacecraft = np.less(distance_to_spacecraft, burst_diamaters/2)
        n[i] = np.sum(close_to_spacecraft)
    total_detected = np.sum(n)
    cdf = np.array([np.sum(n[i:])/total_detected for i in range(len(n))])
    return cdf

def metroplis(start, target, proposal, niter, nburn=0, 
            thin=1, verbose=False, log_likelihood=False):
    """
    This function implements the Metropolis–Hastings sampler.
    Does not use the log-Likelihood yet.
    """
    niter = int(niter)
    current = start
    post = -1E31*np.ones((niter, len(start)), dtype=float)

    # Determine if the target will be comparing the 
    # likelihood or log-likelihood.
    if log_likelihood:
        compare_value = 0
    else:
        compare_value = 1

    for i in progressbar.progressbar(range(niter)):
        proposed = proposal(current)
        p = min(target(proposed)/target(current), compare_value)
        if np.random.random() < p:
            current = proposed
        post[i, :] = current
        if verbose:
            print(f'Current values {current}')
    return post[int(nburn)::thin, :]

def gaus_likelihood(p, x, y, niter):
    """ Gaussian likelihood. """
    C = (np.std(y)*np.sqrt(2*np.pi))
    n_a = int(p[0]*niter)
    burst_diameters = np.concatenate((
        p[1]*np.ones(n_a),
        p[2]*np.ones(niter-n_a)
        ))
    #print(burst_diameters)
    y_model = mc_brute_vectorized(burst_diameters, 
                            bins=x, n_bursts=niter)
    args = sum([(y_i - y_model_i)**2 
                for (y_i, y_model_i) in zip(y, y_model)])
    return np.exp(-0.5*args/LIKELIHOOD_ERROR**2)/C

def proposal(p, proposal_jump=[0.01, 5, 5]):
    """ 
    Generate a new proposal, or "guess" for the MCMC to try next. 
    The new proposed value is picked from a Normal, centered on 
    the old value given by p. The proposal_jump array specifies 
    the standard deviation of the possible jumps from the prior
    parameter value.
    """
    new_vals = np.array([scipy.stats.norm(loc=p_i, scale=jump_i).rvs() 
                    for p_i, jump_i in zip(p, proposal_jump)])  
    # Keep the mixing probability between 0 and 1 inclusive.
    new_vals[new_vals < 0] = 0 # Keep the microburst sizes positive as well.
    #if new_vals[0] < 0: new_vals[0] = 0
    if new_vals[0] > 1: 
        new_vals[0] = 1
    return new_vals

if __name__ == '__main__':
    # Load the CDF data to model
    cdf_data = pd.read_csv(CDF_DATA_PATH)
        
    # Specify priors for the one fixed-sized microburst model.
    # Two parameter model. First parameter is the mixing term, and second and third 
    # parameters are the two microburst sizes.
    # prior = [scipy.stats.halfnorm(loc=0, scale=60)]
    prior = [
            scipy.stats.halfnorm(loc=0, scale=0.3), 
            scipy.stats.norm(loc=100, scale=50),
            scipy.stats.norm(loc=30, scale=20) 
            ]
    # Initial guess on the microburst size.
    start = [prior_i.rvs() for prior_i in prior]

    # The target function. If probability is higher, take the new value given from proposal. Else do the Metroplis thing where you draw a random number between 
    # 0 and 1 and compare to the target value (which will be less than 1).
    niter = 100000
    target = lambda p: gaus_likelihood(p, cdf_data['Separation [km]'], 
                                cdf_data['CDF'], niter)*np.prod(
                                [prior_i.pdf(p_i) for prior_i, p_i in zip(prior, p)])

    if OVERWRITE or (not os.path.exists(SAVE_PATH)):
        trace = metroplis(start, target, proposal, niter, 
                                nburn=10000, thin=1, verbose=False)
        # Save data
        df = pd.DataFrame(data=trace, columns=['a', 'd0', 'd1'])
        df.to_csv(SAVE_PATH, index=False)

    else:
        print('Data already saved. Aborting MCMC.')
        df = pd.read_csv(SAVE_PATH)
        
    print(df.quantile([0.025, 0.5, 0.975]))

    ### PLOTTING CODE ###
    fig = plt.figure(figsize=(10, 8))
    gs = gridspec.GridSpec(3, 3)
    ax = np.zeros((2, 3), dtype=object)
    for row in range(ax.shape[0]):
        for column in range(ax.shape[1]):
            ax[row, column] = plt.subplot(gs[row, column])
    bx = plt.subplot(gs[2, :])

    N = df.shape[0]
    colors = ['g', 'r', 'b']
    ax[0,0].plot(np.arange(N)/1E4, df.a, c='k')
    ax[1,0].hist(df.a, density=True, bins=np.linspace(0, 1), color='k')
    ax[0,0].set(xlabel=r'Iteration x $10^4$', ylabel='a trace')
    ax[1,0].plot(np.linspace(0, 1), prior[0].pdf(np.linspace(0, 1)))
    ax[1,0].set(xlabel='a', ylabel='a posterior PD')

    ax[0,1].plot(np.arange(N)/1E4, df.d0, c='k')
    ax[0,1].set(xlabel=r'Iteration x $10^4$', ylabel=r'$d_0$ trace')
    ax[1,1].hist(df.d0, density=True, bins=np.linspace(0, 200), color='k')
    ax[1,1].plot(np.linspace(0, 200), prior[1].pdf(np.linspace(0, 200)))
    ax[1,1].set(xlabel=r'$d_0$ [km]', ylabel=r'$d_0$ posterior PD')

    ax[0,2].plot(np.arange(N)/1E4, df.d1, c='k')
    ax[0,2].set(xlabel=r'Iteration x $10^4$', ylabel=r'$d_1$ trace')
    ax[1,2].hist(df.d1, density=True, bins=np.linspace(0, 200), color='k')
    ax[1,2].plot(np.linspace(0, 200), prior[2].pdf(np.linspace(0, 200)))
    ax[1,2].set(xlabel=r'$d_1$ [km]', ylabel=r'$d_1$ posterior PD')

    # Pick 1000 traces to analyze further to make plots.
    rand_ind = np.random.choice(np.arange(N), size=1000)
    y_model = np.nan*np.zeros((len(rand_ind), 
                            len(cdf_data['Separation [km]'])))

    # Plot 100 random traces on top of the data.
    N_plot = 100
    j = 0
    for _, row in df.loc[rand_ind, :].iterrows():
        n_a = int(row.a*niter)
        burst_diameters = np.concatenate((
            row[1]*np.ones(n_a),
            row[2]*np.ones(niter-n_a)
            ))
        #print(burst_diameters)
        y_model[j, :] = mc_brute_vectorized(burst_diameters, 
                                bins=cdf_data['Separation [km]'])
        j += 1

    for i in range(N_plot):
        bx.plot(cdf_data['Separation [km]'], y_model[i,:], c='grey', alpha=0.2)
    bx.plot(cdf_data['Separation [km]'], cdf_data['CDF'], c='k')

    # Find the mean and 95% interval of the 1000 curves.
    quartiles = [2.5, 50, 97.5]
    y_quartile = np.percentile(y_model, quartiles, axis=0)

    for i, q in enumerate(y_quartile):
        bx.plot(cdf_data['Separation [km]'], q, c=colors[i], 
                label=f'{quartiles[i]}')

    ax[0, 1].set_title('Two microburst population MCMC model\n'
                       r'$pdf = a \delta(d-d_0) + (1-a) \delta(d-d_1)$')
    bx.set(xlabel='Spacecraft separation [km]', ylabel='F(d)')
    gs.tight_layout(fig)
    plt.show()
