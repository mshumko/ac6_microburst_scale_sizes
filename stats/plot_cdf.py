import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser

converters = {0:dateutil.parser.parse, 
            -1:dateutil.parser.parse, 
            -2:dateutil.parser.parse}

bins = np.arange(0, 70, 5)
frac = np.nan*np.zeros(len(bins)-1)
CC_thresh = 0.8

dtypes = (object, float, float, float, float, float, float, float, 
        float, float, float, float, float, float, float, object, object)
data = np.genfromtxt('coincident_microburst_test.csv', delimiter=',',
        names=True, dtype=dtypes)#, converters=converters)
# years = np.array([t.year for t in data['dateTime']])
# months = np.array([t.month for t in data['dateTime']])

for i, (lower_edge, upper_edge) in enumerate(zip(bins[:-1], bins[1:])):
        # Find microbursts in bin
        idsep = np.where(
                        (data['Dist_Total'] > lower_edge) & 
                        (data['Dist_Total'] <= upper_edge) &
                        #np.logical_not((years == 2015) & (months == 4)) &
                        #(np.abs(data['Lm_OPQ']) < 8) & (data['dos1rate'] > 100)
                        (data['time_cc'] > data['space_cc'])
                        )[0]
        idCoincident = np.where((data['time_cc'][idsep] >= CC_thresh))[0]
        if len(idsep):
                frac[i] = len(idCoincident)/len(idsep)

        if lower_edge >= 85 and lower_edge <= 90:
                print(lower_edge, upper_edge, data['dateTime'][idsep], '\n')

pdf = np.convolve([-1, 1], frac, mode='valid')

fig, ax = plt.subplots(1, sharex=True)
ax.bar(np.convolve([0.5, 0.5], bins, mode='valid'), frac, width=(bins[1]-bins[0])*0.8)
ax.set_ylabel('CDF')

# ax[1].bar(bins[1:-1], pdf, width=(bins[1]-bins[0])*0.8)
# ax[1].set_ylabel('PDF')

plt.show()
