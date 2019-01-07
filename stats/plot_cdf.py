import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser

# converters = {0:dateutil.parser.parse, 
#             -1:dateutil.parser.parse, 
#             -2:dateutil.parser.parse}

bins = np.arange(0, 50, 5)
frac = np.nan*np.zeros(len(bins)-1)
CC_thresh = 0.7

dtypes = (object, float, float, float, float, float, float, float, 
        float, float, float, float, float, float, float, object, object)
data = np.genfromtxt('coincident_microburst_test.csv', delimiter=',',
                    names=True, dtype=dtypes)

for i, (lower_edge, upper_edge) in enumerate(zip(bins[:-1], bins[1:])):
    # Find microbursts in bin
    idsep = np.where((data['Dist_Total'] > lower_edge) & 
            (data['Dist_Total'] <= upper_edge))[0]
    idCoincident = np.where((data['time_cc'][idsep] >= CC_thresh) & 
                    (data['time_cc'][idsep] > data['space_cc'][idsep]))[0]
    frac[i] = len(idCoincident)/len(idsep)

pdf = np.convolve([-1, 1], frac, mode='valid')

fig, ax = plt.subplots(2, sharex=True)
ax[0].bar(bins[:-1], frac, width=(bins[1]-bins[0])*0.8)
ax[0].set_ylabel('CDF')

ax[1].bar(bins[1:-1], pdf, width=(bins[1]-bins[0])*0.8)
ax[1].set_ylabel('PDF')

plt.show()