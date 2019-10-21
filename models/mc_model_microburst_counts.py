# This program contains the microburst Monte Carlo model that instead of 
# calculating the modeled microburst CDF, calculates the number of 
# microbursts observed in each separation bin given AC6's uneven 
# sampling in separation.

import numpy as np
import pandas as pd

def mc_vectorized(burst_diamaters, weights, n_obs,
                n_bursts=int(1E5), grid_size_km=200,
                s_bins_km=np.arange(0, 100, 5)):
    """ 
    This function tallies events that were observed by both spacercaft in each
    AC6 separation bin. The number of microbursts observed in each separation 
    bin is weighted by the relative fraction of simultanous 10Hz samples taken
    for each separation bin.
    """
    n_model = np.zeros(len(s_bins_km))
    
    # If supplied a single-valued diameter burst_diamaters. If supplied an array, 
    # I implicityly assume the diameters are distributed according to 
    # some PDF distribution.
    if not hasattr(burst_diamaters, '__len__'):
        burst_diamaters = burst_diamaters*np.ones(n_bursts)
    # Randomly generate n_burst microburst centers.
    burst_centers = np.random.uniform(-grid_size_km, grid_size_km, size=(n_bursts, 2))
    
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
    for i, d in enumerate(s_bins_km):
        spacecraft_location = np.zeros_like(burst_centers)
        spacecraft_location[:, 1] = d
        distance_to_spacecraft = np.linalg.norm(burst_centers-spacecraft_location, 
                                            axis=1)
        # Where close_to_origin is True, the microburst intersected the origin.
        # We only need to loop over the True values.
        close_to_spacecraft = np.less(distance_to_spacecraft, burst_diamaters/2)
        n_model[i] = np.sum(close_to_spacecraft)
    total_detected = np.sum(n_model)
    return n_model*(n_obs/total_detected)