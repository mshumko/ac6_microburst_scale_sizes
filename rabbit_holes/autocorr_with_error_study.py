# Test code to see how the range 
import numpy as np
import matplotlib.pyplot as plt

N = 1000
y = np.array([250.001, 270.001, 270.001, 320.002, 370.002, 330.002, 329.969,
       310.002, 270.001])
cc = np.zeros(N)

for i in range(N):
    y_noise = np.random.poisson(y)
    cc_unnormalized = np.correlate(y-y.mean(), y_noise-y_noise.mean(), mode='valid')
    cc[i] = cc_unnormalized/np.sqrt(len(y)*len(y)*np.var(y_noise)*np.var(y))

fig, ax = plt.subplots()
ax.hist(cc, range=[.75, 1])
ax.text(0.01, 0.99, 'N={}'.format(N), transform=ax.transAxes, va='top')
ax.set_ylabel('Count Rate'); ax.set_xlabel('Cross-correlation')
plt.show()