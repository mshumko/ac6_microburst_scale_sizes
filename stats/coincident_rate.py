import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
import datetime
import dateutil.parser
import csv
import itertools
import sys
import os 

sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

class OccuranceRate:
    def __init__(self, sc_id, date, catV, catPath=None):
        """
        For a given spacecraft on a given day, this class calculates
        the microburst occurance rate (with an assumption about width.)
        """
        self.sc_id = sc_id
        self.date = date
        self._load_sc_data()

        if catPath is None:
            catPath = ('/home/mike/research/ac6-microburst-scale-sizes'
                      '/data/microburst_catalogues')
        self._load_catalog(catPath, catV)
        return

    def radBeltIntervals(self, lbound=4, ubound=8):
        """
        This method finds all of the radiation belt passes in a day 
        and separates it into start/end indicies.
        """
        # Find the start/end indicies of each rad belt pass.
        L = np.abs(self.data['Lm_OPQ'])
        idL = np.where((L > lbound) & (L < ubound))[0]
        conv = np.convolve([1, -1], idL, mode = 'valid') - 1
        consecutiveFlag = np.where(conv != 0)[0] + 1
        startInd = np.insert(consecutiveFlag, 0, 0)
        endInd = np.insert(consecutiveFlag, 
                len(consecutiveFlag), len(idL)-1)
        # Save pass indicies to self.intervals.
        self.intervals = np.zeros((len(startInd), 2), dtype=int)
        for i, (i_s, i_e) in enumerate(zip(startInd, endInd)):
            self.intervals[i, :] = [idL[i_s], idL[i_e-1]]
        return 

    def occurance_rate(self, mode='static', ):
        """ 
        This method calculates the microburst occurance rates
        for each pass assuming a microburst width assumption.
        The width is either a fixed value in seconds or 
        another class method.
        """
        if mode == 'static':
            width = self.static_width(0.5)
        nDetTime = date2num(self.cat['dateTime'])
        self.rates = -1*np.ones(self.intervals.shape[0], 
                                dtype=float)
        for i, (i_start, i_end) in enumerate(self.intervals):
            # Find microbursts in ith pass
            startTime = date2num(self.data['dateTime'][i_start])
            endTime = date2num(self.data['dateTime'][i_end])
            idt = np.where((nDetTime > startTime) & 
                            (nDetTime < endTime))[0]
            #print('Pass', i, 'len(idt)=', len(idt))
            pass_duration = (self.data['dateTime'][i_end] -
                self.data['dateTime'][i_start]).total_seconds()
            self.rates[i] = 100*width*len(idt)/pass_duration  
        return
        
    def static_width(self, width):
        """ Returns a constant width """
        return width
        
    def prominence_width(self):
        raise NotImplementedError('Topological prominence '
                                    'not working yet')
        return

    def _load_sc_data(self):
        """ Loads AC-6 10 Hz data """
        print('Loading AC6-{} data from {}'.format(
            self.sc_id.upper(), self.date.date()))
        self.data = read_ac_data.read_ac_data_wrapper(
                                self.sc_id, self.date)
        return

    def _load_catalog(self, catPath, v):
        """ 
        Loads the microburst catalog of version v spacecraft given 
        by self.sc_id. 
        """
        fPath = os.path.join(catPath, 
            'AC6{}_microbursts_v{}.txt'.format(self.sc_id.upper(), v))
        print('Loading catalog,', 'AC6{}_microbursts_v{}.txt'.format(
            self.sc_id.upper(), v))
        with open(fPath) as f:
            r = csv.reader(f)
            keys = next(r)
            rawData = np.array(list(r))
        self.cat = {}
        for (i, key) in enumerate(keys):
            self.cat[key] = rawData[:, i]

        # Now convert the times array(s) to datetime, 
        # and all others except 'burstType' to a float.
        timeKeys = itertools.filterfalse(
            lambda x:'dateTime' in x or 'burstType' in x, self.cat.keys())
        for key in timeKeys:
                self.cat[key] = self.cat[key].astype(float)
        for key in filter(lambda x: 'dateTime' in x, self.cat.keys()):
            self.cat[key] = np.array([dateutil.parser.parse(i) 
                for i in self.cat[key]])
        return

if __name__ == '__main__':
    o = OccuranceRate('A', datetime.datetime(2016, 10, 14), 3)
    o.radBeltIntervals()
    o.occurance_rate()
    startTimes = o.data['dateTime'][o.intervals[:, 0]]
    _, ax = plt.subplots(2, sharex=True)
    ax[0].scatter(startTimes, o.rates)
#    plt.bar(startTimes, o.rates, width=passDur align='left')
    #ax[0].set_ylabel('Occurance rate\n(microbursts/pass)')
    ax[0].set_ylabel('Microbust duty cycle\n(microburst/pass time')
    
    idt = np.where((o.cat['dateTime'] > startTimes[0]) & 
                    (o.cat['dateTime'] < startTimes[-1]))[0]
    ax[1].scatter(o.cat['dateTime'][idt], o.cat['AE'][idt])
    ax[1].set_ylabel('AE (nT)')
    ax[1].set_xlabel('UTC')
    plt.show()
