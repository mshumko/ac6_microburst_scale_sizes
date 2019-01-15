import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import dateutil.parser
import csv

openPath = 'coincident_microburst_test.csv'
savePath = 'coincident_microburst_test_v2.csv'
dtype = [datetime] + 14*[float] + 2*[datetime]
converters = {0:dateutil.parser.parse,
            -2:dateutil.parser.parse,
            -1:dateutil.parser.parse}

data = np.genfromtxt(openPath, delimiter=',', names=True, dtype=dtype, converters=converters)
# Sort data
data_sorted = sorted(data, key=lambda x:x[0])
# Remove duplicates
t = np.array([row[0] for row in data_sorted])
dt = np.array([i.total_seconds() for i in (t[1:] - t[:-1])])

# Save file
idt = np.where(dt > 0.2)[0] # Non-duplicate values
with open(savePath, 'w') as f:
    w = csv.writer(f, delimiter=',')
    w.writerow(data.dtype.names)
    for i in idt:
        w.writerow(data_sorted[i])

