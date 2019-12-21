import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

CATALOG_DIR = ('/home/mike/research/ac6_microburst_scale_sizes/'
                'data/coincident_microbursts_catalogues')
CATALOG_NAME = 'AC6_coincident_microbursts_2scc_v9_sorted.txt'

df = pd.read_csv(os.path.join(CATALOG_DIR, CATALOG_NAME))
df = df.dropna(subset=['peak_width_A', 'peak_width_B'])

fig, ax = plt.subplots(1, 2, sharex=True)

ax[0].hist(df.peak_width_A)
ax[1].hist(df.peak_width_B)

n_A = df.peak_width_A.shape[0]
n_B = df.peak_width_B.shape[0]

n_A_gt1 = df[df.peak_width_A >= 1].shape[0]
n_B_gt1 = df[df.peak_width_B >= 1].shape[0]

print(n_A, n_A_gt1, n_B, n_B_gt1)

ax[0].set(ylabel='Number of microbursts', xlabel='Peak width [s]', title='AC6A')
ax[1].set(xlabel='Peak width [s]', title='AC6B')

plt.show()