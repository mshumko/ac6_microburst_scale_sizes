import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser

converters = {0:dateutil.parser.parse, 
            -1:dateutil.parser.parse, 
            -2:dateutil.parser.parse}

bins = np.arange(0, 200, 5)
frac = np.nan*np.zeros(len(bins)-1)
num = np.nan*np.zeros(len(bins)-1)
CC_thresh = 0.8

dtypes = (object, float, float, float, float, float, float, float, 
        float, float, float, float, float, float, float, object, object)
data = np.genfromtxt('coincident_microburst_test_v2.csv', delimiter=',',
        names=True, dtype=dtypes)#, converters=converters)
# years = np.array([t.year for t in data['dateTime']])
# months = np.array([t.month for t in data['dateTime']])

curtain_thresh = 0.1

fig, ax = plt.subplots(3, figsize=(8, 10))

for i, (lower_edge, upper_edge) in enumerate(zip(bins[:-1], bins[1:])):
        # Find microbursts in bin
        idsep = np.where(
                        (data['Dist_Total'] > lower_edge) & 
                        (data['Dist_Total'] <= upper_edge) &
                        # Filter by L shell
                        (np.abs(data['Lm_OPQ']) < 8) &
                        (np.abs(data['Lm_OPQ']) > 4) &
                        (data['time_cc'] > data['space_cc']+curtain_thresh) &
                        # Geographic filter to pass data outside of the US.
                        (((data['lon'] > -60) | (data['lon'] < -140)) |
                        ((data['lat'] > 60) | (data['lat'] < -60)))
                        )[0]
        idCoincident = np.where((data['time_cc'][idsep] >= CC_thresh))[0]
        if len(idsep):
                print('lower_edge =', lower_edge, 'upper_edge =', upper_edge, 'len(idsep) =', len(idsep))
                frac[i] = len(idCoincident)/len(idsep)
                num[i] = len(idsep)
        ax[2].scatter(data['lon'][idsep], data['lat'][idsep])

        # if lower_edge >= 85 and lower_edge <= 90:
        #         print(lower_edge, upper_edge, data['dateTime'][idsep], '\n')

pdf = np.convolve([-1, 1], frac, mode='valid')

#fig, ax = plt.subplots(2)
ax[0].bar(np.convolve([0.5, 0.5], bins, mode='valid'), frac, width=(bins[1]-bins[0])*0.8)
ax[0].set_ylabel('CDF')

ax[1].bar(np.convolve([0.5, 0.5], bins, mode='valid'), num, width=(bins[1]-bins[0])*0.8)
ax[1].set_ylabel('Number of detections')

# ax[1].bar(bins[1:-1], pdf, width=(bins[1]-bins[0])*0.8)
# ax[1].set_ylabel('PDF')

#ax[1].scatter(data['lon'][idsep], data['lat'][idsep])

plt.show()
