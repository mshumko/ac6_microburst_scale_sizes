# This script calls the 
import matplotlib.pyplot as plt

import sys
sys.path.append('/home/mike/research/ac6-microburst-scale-sizes/flash-curtain-determination')
import replace_error_sep_lags

r = replace_error_sep_lags.ReplaceErrorVals('/home/mike/research/ac6/AC6_Separation.csv')
r.loadSeprationFile()

# Plot the lifetime AC-6 separation
plt.plot(r.sepDict['Date/Time'], r.sepDict['In-Track Separation [km]'])
plt.xlabel('Date')
plt.ylabel('Separation [km]')
plt.title('AC-6 In-Track Separation')
plt.tight_layout()
plt.show()