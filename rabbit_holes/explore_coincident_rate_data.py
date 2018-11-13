# This script explores the coincident rate data.
import numpy as np
import matplotlib.pyplot as plt
import os

data_name = 'microburst_catalog.csv'
directory = ('/home/mike/research/ac6_microburst_scale_sizes/'
            'data/coincident_microbursts_catalogues')
data_path = os.path.join(directory, data_name)

# Load data
d = np.genfromtxt(data_path, delimiter=',', names=True, usecols=range(1,5))
dist = np.array([i[0] for i in d])
ratio = np.array([i[-1] for i in d])

# Replace error values
# ratio = np.nan_to_num(ratio)
ratio[np.where(np.isinf(ratio))[0]] = -1
ratio[np.where(np.isnan(ratio))[0]] = -1

# Scatter plot the results
plt.scatter(dist, ratio)
plt.xlim(0, 150); plt.ylim(0, 500)
plt.show()