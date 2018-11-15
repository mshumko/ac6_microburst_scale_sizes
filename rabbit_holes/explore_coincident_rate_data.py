# This script explores the coincident rate data.
import numpy as np
import matplotlib.pyplot as plt
import os

data_name = 'microburst_catalog_v2.csv'
directory = ('/home/mike/research/ac6_microburst_scale_sizes/'
            'data/coincident_microbursts_catalogues')
data_path = os.path.join(directory, data_name)

# Load data
d = np.genfromtxt(data_path, delimiter=',', names=True, usecols=range(1,5))
dist = np.array([i[0] for i in d])
microburst_occurance = np.array([i[1]-i[2] for i in d])

# Replace error values
# ratio = np.nan_to_num(ratio)
#ratio[np.where(np.isinf(ratio))[0]] = np.nan
#ratio[np.where(np.isnan(ratio))[0]] = -1

bins = np.arange(0, 151, 10)
hist = {i:microburst_occurance[np.where((dist > i) & (dist <= j) & (~np.isnan(microburst_occurance)
& (microburst_occurance > 0)))[0]] for (i, j) in 
        zip(bins[:-1], bins[1:]) }

labels, data = [*zip(*hist.items())]
plt.boxplot(data)
plt.xticks(range(1, len(labels) + 1), labels)

# Scatter plot the results
#plt.scatter(dist, ratio)
#plt.xlim(0, 150); plt.ylim(0, 500)
#plt.title('AC6 coincident/chance microburst ratio\ncatalog version =1 | CC_thresh = 0.7')
plt.xlabel('Total_separation [km]')
#plt.ylabel('coincident/chance microburst ratio')
plt.show()
