import numpy as np
import matplotlib.pyplot as plt
import dateutil.parser

write_times_to_file = False
write_cdf_to_file = False

version = 4
catPath = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
        'coincident_microbursts_catalogues/'
       'AC6_coincident_microbursts_v{}.txt'.format(version))

converters = {0:lambda t: dateutil.parser.parse(t.decode()), 
            -1:lambda t: dateutil.parser.parse(t.decode()), 
            -2:lambda t: dateutil.parser.parse(t.decode())}

bins = np.arange(2, 100, 3)
frac = np.nan*np.zeros(len(bins)-1)
num = np.nan*np.zeros(len(bins)-1)
CC_thresh = 0.8

dtypes = [object] + [float]*15 + 2*[object]
data = np.genfromtxt(catPath, delimiter=',',
        names=True, dtype=dtypes)#, converters=converters)


curtain_thresh = 0

fig, ax = plt.subplots(3, figsize=(6, 10))

for i, (lower_edge, upper_edge) in enumerate(zip(bins[:-1], bins[1:])):
        # Find microbursts in bin
        idsep = np.where(
                        #(data['Dist_Total'] > lower_edge) 
                        #  & 
                        (data['Dist_Total'] <= upper_edge) 
                        & # Filter by significance above the 10% baseline.
                        (data['peak_std'] > 2)
                        &
                        # Filter by L shell
                        (np.abs(data['Lm_OPQ']) < 8) 
                        &
                        (np.abs(data['Lm_OPQ']) > 4) 
                        &
                        # Filter by a temporal CC > spatial CC (and threshold).
                        (data['time_cc'] > data['space_cc']+curtain_thresh) 
                        &
                        # Geographic filter to filter data inside the US.
                        (
                        ((data['lon'] > -60) | (data['lon'] < -140)) |
                        ((data['lat'] > 70) | (data['lat'] < 15))
                        ) 
                        &
                        # Geographic filter to filter data inside the SAA.
                        (
                        ((data['lon'] > 30)  | (data['lon'] < -116)) |
                        ((data['lat'] < -90) | (data['lat'] > 0))
                        )
                        )[0]
        # Find the subset of idsep events where 
        # temporal CC > threshold.                
        idCoincident = np.where((data['time_cc'][idsep] >= CC_thresh))[0]
        
        if len(idsep):
            print('lower_edge =', lower_edge, 'upper_edge =', upper_edge, 'len(idsep) =', len(idsep))
            frac[i] = len(idCoincident)/len(idsep)
            num[i] = len(idsep)
        ax[2].scatter(data['lon'][idsep], data['lat'][idsep], color='k')

        if write_times_to_file and lower_edge == 80 and upper_edge == 85:
            with open('test_times.csv', 'w') as f:
                for ti in data['dateTime'][idsep]:
                    f.write(str(ti) + '\n')
            #print(lower_edge, upper_edge, data['dateTime'][idsep], '\n')

        if write_cdf_to_file:
            saveArr = np.stack([bins[:-1], frac])
            np.savetxt('cdf_vs_sep.csv', saveArr.T, fmt='%.3e', delimiter=', ', 
                        header='sep_km, cdf')

#pdf = np.convolve([-1, 1], frac, mode='valid')
frac /= max(frac)
ax[0].bar(np.convolve([0.5, 0.5], bins, mode='valid'), frac, width=(bins[1]-bins[0])*0.8, color='k')
ax[0].set_ylabel('CDF (Normalized to 1)'); #ax[0].set_ylim(bottom=0.15)
ax[0].set_xlabel('AC-6 total separation [km]')

# ax[0].bar(np.convolve([0.5, 0.5], bins[:-2], mode='valid'), 
#         np.convolve([-1, 1], frac[1:], mode='valid'), width=(bins[1]-bins[0])*0.8)
# ax[0].set_ylabel('PDF'); #ax[0].set_ylim(bottom=0.15)

ax[1].bar(np.convolve([0.5, 0.5], bins, mode='valid'), num, width=(bins[1]-bins[0])*0.8, color='k')
ax[1].set_ylabel('cumulative number of detections\ni.e. number of detections in bin 0 to x km')
ax[1].set_xlabel('AC-6 total separation [km]')

ax[-1].set(xlabel='Lon', ylabel='Lat')
plt.tight_layout()
plt.show()
