# Test code to see how the range 
import numpy as np
import matplotlib.pyplot as plt

N = 1000
y = np.array([250.001, 270.001, 270.001, 320.002, 370.002, 330.002, 329.969,
       310.002, 270.001])/10
cc = np.zeros(N)

error = ''

for i in range(N):
    y_noise = np.random.poisson(y)

    # If error mode is 'full', then assume an error from the count rates and 
    # mean subtraction. Add these uncertanties in quadreture.
    if error == 'full':
        residuals = y - y_noise
        y_noise = y + np.sign(residuals)*np.sqrt(residuals**2 + np.var(y))
        
    cc_unnormalized = np.correlate(y-y.mean(), y_noise-y_noise.mean(), mode='valid')
    cc[i] = cc_unnormalized/np.sqrt(len(y)*len(y)*np.var(y_noise)*np.var(y))

fig, ax = plt.subplots()
ax.hist(cc, range=[0, 1])
ax.text(0.01, 0.99, 'N={}'.format(N), transform=ax.transAxes, va='top')
ax.set_ylabel('Count Rate'); ax.set_xlabel('Cross-correlation')
plt.show()
