# This program run a MCMC simulation of a fixed-sized microburst population

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import os

import pandas as pd 
import progressbar

GRID_SIZE = 200
PRIOR = 'uniform'
# csv save data path. Will NOT overwrite if it already exists!
SAVE_PATH = ('/home/mike/research/ac6_microburst_scale_sizes/models/mcmc_traces'
            f'/mcmc_one_size_{PRIOR}_trace.csv')
CDF_DATA_PATH = ('/home/mike/research/ac6_microburst_scale_sizes'
            '/data/microburst_cdf_pdf_norm_v3.csv')
PAPER_PLOT = True
OVERWRITE = False
    
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
    y_model = mc_brute_vectorized(p[0], bins=x)
    args = sum([(y_i - y_model_i)**2 
                for (y_i, y_model_i) in zip(y, y_model)])
    return np.exp(-0.5*args/np.var(y))/C

def proposal(p, proposal_jump=[5]):
    """ 
    Generate a new proposal, or "guess" for the MCMC to try next. 
    The new proposed value is picked from a Normal, centered on 
    the old value given by p. The proposal_jump array specifies 
    the standard deviation of the possible jumps from the prior
    parameter value.
    """
    new_vals = np.array([scipy.stats.norm(loc=p_i, scale=jump_i).rvs() 
                    for p_i, jump_i in zip(p, proposal_jump)])         
    return new_vals

if __name__ == '__main__':
    # Load the CDF data to model
    cdf_data = pd.read_csv(CDF_DATA_PATH)
        
    # Specify prior for the one fixed-sized microburst model.
    if PRIOR == 'uniform':
        prior = [scipy.stats.uniform(0, 200)]
    else:
        prior = [scipy.stats.halfnorm(loc=0, scale=60)]
    # Initial guess on the microburst size.
    start = [prior_i.rvs() for prior_i in prior]

    # The target function. If probability is higher, take the new value given from proposal. Else do the Metroplis thing where you draw a random number between 
    # 0 and 1 and compare to the target value (which will be less than 1).
    target = lambda p: Likelihood(p, cdf_data['Separation [km]'], 
                                cdf_data['CDF'])*np.prod(
                                [prior_i.pdf(p_i) for prior_i, p_i in zip(prior, p)])
    niter = 100000

    if OVERWRITE or (not os.path.exists(SAVE_PATH)):
        trace = metroplis(start, target, proposal, niter, 
                                nburn=10000, thin=1, verbose=False)
        # Save data
        df = pd.DataFrame(data=trace, columns=['d'])
        df.to_csv(SAVE_PATH, index=False)
    else:
        print('Data already saved. Aborting MCMC.')
        df = pd.read_csv(SAVE_PATH)
    
    print(df.quantile([0.025, 0.5, 0.975]))

    ### PLOTTING CODE ###
    if not PAPER_PLOT:
        colors = ['r', 'g', 'b']
        labels = ['2.5%', '50%', "97.5%"]
        _, ax = plt.subplots(3, 1, figsize=(8, 9))
        ax[0].plot(df.d, c='k')
        ax[1].hist(df.d, density=True, bins=np.arange(0, 200), color='k', label='posterior')
        ax[1].plot(np.linspace(0, 200), prior[0].pdf(np.linspace(0, 200)), c='c', label='prior')
        ax[2].plot(cdf_data['Separation [km]'], cdf_data['CDF'], c='k', label='CDF data')
        for i, size in enumerate(np.percentile(df.d, [2.5, 50, 97.5])):
            ax[2].plot(cdf_data['Separation [km]'], 
                        mc_brute_vectorized(size, bins=cdf_data['Separation [km]']), 
                        c=colors[i])
            ax[1].axvline(size, c=colors[i], label=labels[i])

        ax[0].set_title('One size microburst MCMC model | prior ~ U(0, 200)')
        ax[0].set(xlabel='Iteration', ylabel='Microburst diameter trace'); 

        ax[1].set(xlabel='microburst diameter [km]', ylabel='Probability probability'); 
        ax[1].legend()
        ax[2].legend()
        ax[2].set(xlabel='microburst diameter [km]', ylabel='Microburst fraction')

        plt.tight_layout()
        plt.show()
    else:
        colors = ['r', 'g', 'b']
        labels = ['2.5 %', '50 %', "97.5 %"]
        _, ax = plt.subplots(2, 1, figsize=(6, 5))
        ax[0].hist(df.d, density=True, bins=np.arange(0, 200), color='k', label='_nolegend_')
        ax[0].plot(np.linspace(0, 200), prior[0].pdf(np.linspace(0, 200)), c='c')

        # Plot 100 random traces on top of the data.
        N_plot = 100
        rand_ind = np.random.choice(np.arange(df.shape[0]), size=N_plot)
        y_model = np.nan*np.zeros((len(rand_ind), 
                                len(cdf_data['Separation [km]'])))

        for _, row in df.loc[rand_ind, :].iterrows():
            burst_diameters = row.d
            y_model = mc_brute_vectorized(burst_diameters, 
                                    bins=cdf_data['Separation [km]'])
            ax[1].plot(cdf_data['Separation [km]'], y_model, c='grey', alpha=0.3)

        # plot the AC6 data
        ax[1].plot(cdf_data['Separation [km]'], cdf_data['CDF'], c='k', label='AC6 F(s)')
        # Plot the quantiles
        for i, size in enumerate(np.percentile(df.d, [2.5, 50, 97.5])):
            ax[1].plot(cdf_data['Separation [km]'], 
                        mc_brute_vectorized(size, bins=cdf_data['Separation [km]']), 
                        c=colors[i], label=labels[i])
            ax[0].axvline(size, c=colors[i])

        ax[0].set_title('One microburst size MCMC model\n'
                        r'$pdf = \delta(s-d)$')

        ax[0].set(xlabel='d [km]', ylabel='posterior PD', xlim=(26, 140)); 
        ax[1].legend()
        ax[1].set(xlabel='AC6 separation (s) [km]', ylabel='F(s)', xlim=(0, 90))

        plt.tight_layout()
        plt.show()
