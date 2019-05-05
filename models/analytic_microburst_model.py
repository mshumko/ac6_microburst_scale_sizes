# Now that Paul motivated me to do these integrals.
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

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


if __name__ == '__main__':
    # New method
    diameter_km = 40 
    s_array = np.arange(0, 50)
    cdf = one_size_cdf(diameter_km, s_array)

    # Old method
    A_array = A(diameter_km/2, s_array)

    cdf_2 = np.array([np.nansum(A_array[i:])/np.nansum(A_array) for i in range(len(s_array))])

    plt.plot(s_array, cdf)
    plt.plot(s_array, cdf_2, '--')
    plt.show()