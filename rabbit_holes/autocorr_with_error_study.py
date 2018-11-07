# Test code to see how the range 
import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser

N = 1000
# y = np.array([250.001, 270.001, 270.001, 320.002, 370.002, 330.002, 329.969,
#        310.002, 270.001])/10

dataPath = ('/home/mike/research/ac6_microburst_scale_sizes/'
            'data/train/ac6a_training_data.csv')
data = np.genfromtxt(dataPath, delimiter=',', usecols=np.arange(1, 19))
cc = np.zeros(10)
cc_bins = np.arange(0, 1.1, 0.1)

error = 'full'

width = 0.5
idx = np.arange(
            len(data[0])//2-10*width//2,
            len(data[0])//2+10*width//2-1,
            dtype=int
                )

for y in data[:, idx]:
    cc_temp = np.zeros(N)
    for i in range(N):
        y_noise = np.random.poisson(y)

        # If error mode is 'full', then assume an error from the count rates and 
        # mean subtraction. Add these uncertanties in quadreture.
        if error == 'full':
            residuals = y - y_noise
            y_noise = y + np.sign(residuals)*np.sqrt(residuals**2 + np.var(y))
            
        cc_unnormalized = np.correlate(y-y.mean(), y_noise-y_noise.mean(), mode='valid')
        cc_temp[i] = cc_unnormalized/np.sqrt(len(y)*len(y)*np.var(y_noise)*np.var(y))
        #cc[i] += cc_unnormalized/np.sqrt(len(y)*len(y)*np.var(y_noise)*np.var(y))
    #print(np.histogram(cc_temp, bins=cc_bins)[0].shape)
    cc += np.histogram(cc_temp, bins=cc_bins)[0]


fig, ax = plt.subplots()
ax.bar(cc_bins[1:], cc, width=0.8*(cc_bins[1] - cc_bins[0]))
ax.text(0.01, 0.99, 'N = {}\nerror_calc = {}\ncc_width= {} s'.format(
        N, error, width), 
        transform=ax.transAxes, va='top')
ax.set_ylabel('Count Rate'); ax.set_xlabel('Cross-correlation')
ax.set_yscale('log'); 
plt.show()
