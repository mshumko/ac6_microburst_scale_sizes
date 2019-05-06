# Now that Paul motivated me to do these integrals.
import numpy as np
import scipy.integrate
import scipy.stats
import matplotlib.pyplot as plt
import time

try:
    import mcmc_maxwell
    COMPARE_MODEL = True
except:
    print('No mcmc_maxwell module found. Will not do a model comparison'
        ' with the MC model')
    COMPARE_MODEL = False

def A(r, s):
    """ 
    Compute the circle-circle intersection area between two circles
    of radius r with circle centers separated by s.
    """
    a = 2*r**2*np.arccos(s/(2*r))
    b = s/2*np.sqrt(4*r**2-s**2)
    return a-b

def one_size_cdf(diameter, s):
    """ 
    Returns the cdf curve for microbursts of radius r,
    evaluated at points s.
    """
    cdf = np.zeros_like(s, dtype=object)
    r = diameter/2
    for i, s_i in enumerate(s):
        result = scipy.integrate.quad(lambda s_i:A(r, s_i) if ~np.isnan(A(r, s_i)) else 0, s_i, np.inf)
        cdf[i] = result[0]
    cdf /= np.max(cdf)
    return cdf

def continuous_cdf(s, dist, max_scale=1000):
    """
    Calculate the microburst CDF assuming a continuous microburst PDF.
    s is the AC6 separation array, dist is the scipy.stats distribution
    evaluated with the distribution parameters, and max_scale is the 
    upper integral bound for microburst sizes. This bound should technically
    be np.inf, but to speed it up choose a value where the microburst PDF
    is really small. 
    """
    cdf = np.zeros_like(s, dtype=object)
    microburst_pdf = lambda x: dist.pdf(x)
    #r = diameter/2
    f = lambda r, s_i:A(r, s_i)*microburst_pdf(r) if ~np.isnan(A(r, s_i)) else 0

    for i, s_i in enumerate(s):
        result = scipy.integrate.dblquad(f, s_i, np.inf, lambda x:0, lambda x:max_scale)
        cdf[i] = result[0]
    cdf /= np.max(cdf)
    return cdf


if __name__ == '__main__':
    # New method
    diameter_km = 40 
    s_array = np.linspace(0, 60, num=20)
    cdf = one_size_cdf(diameter_km, s_array)

    # Old method
    A_array = A(diameter_km/2, s_array)
    cdf_2 = np.array([np.nansum(A_array[i:])/np.nansum(A_array) for i in range(len(s_array))])

    plt.plot(s_array, cdf)
    plt.plot(s_array, cdf_2, '--')
    plt.show()

    if COMPARE_MODEL:
        dist = scipy.stats.maxwell(loc=0, scale=10)
        # Generalized MC model
        n_bursts = int(1E6)
        burst_diameters = dist.rvs(size=n_bursts)
        cdf_mc = mcmc_maxwell.mc_brute_vectorized(2*burst_diameters, 
                            bins=s_array, n_bursts=n_bursts)
        # Fully generalized method
        cdf_analytic = continuous_cdf(s_array, dist)
        #print(cdf_3)

        plt.plot(s_array, cdf_mc, label='MC')
        plt.plot(s_array, cdf_analytic, '--', label='Analytic')
        plt.legend()
        plt.show()

    