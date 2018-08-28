import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates
import sys
import os
import itertools
import dateutil.parser
import csv
import scipy.signal
from datetime import datetime, timedelta

#sys.path.append('/home/mike/research/ac6-microburst-scale-sizes/stats/')
sys.path.append('/home/mike/research/mission-tools/ac6/')
import read_ac_data
#import leo_scale_sizes

class ValidateDetections:
    def __init__(self, sc_id, cPath, saveDir=None):
        """
        This class contains code to validate the microburst detector algorithm
        by plotting the detections made from each spacecraft side by side. No
        attempt is made to match coincident events. 

        For a user-suplied catalog file (path specified by cPath), the 
        catalog file is first read. Then 
        """
        
        # Load catalog
        self.cData = self._load_catalog_(cPath)
        self.sc_id = sc_id
        # Set save directory
        if saveDir is not None:
            self.saveDir = saveDir
        else:
            self.saveDir = ('/home/mike/research/ac6-microburst-scale-sizes/'
                    'plots/validation/{}/'.format(datetime.now().date()))
        return
        
    def plotLoop(self, pltWidth=5):
        """
        This method loops over the detections in the catalog file and plots 
        the dos1 data.
        """
        #currentDate = datetime.date().min
        self._findDateBounds() # Get index bounds for each date.
        
        fig, self.ax = plt.subplots()
        #self.bx = self.ax.twinx()
        
        # Loop over dates.
        for d in self.dateArr:
            print('Making validation plots from:', d[0])
            # Open file
            d[0] = datetime.combine(d[0], datetime.min.time())
            try:
                self.dataA = read_ac_data.read_ac_data_wrapper('A', d[0])
                self.dataB = read_ac_data.read_ac_data_wrapper('B', d[0])
            except AssertionError as err:
                if (('None or > 1 AC6 files found in' in str(err)) or 
                    'File is empty!' in str(err)):
                    continue
                else:
                    raise

            # Loop over detections on d'th day.
            #zCenterTimes = self.cData['dateTime'][d[1]:d[2]]
            for i, tA in enumerate(self.cData['dateTime'][d[1]:d[2]]): 
                # Plot dos1rate
                pltRange = [tA-timedelta(seconds=pltWidth), 
                            tA+timedelta(seconds=pltWidth)]
                flag = self.plot(self.dataA['dateTime'], 
                    self.dataA['dos1rate'], self.dataB['dateTime'], 
                    self.dataB['dos1rate'], pltRange)
                if flag == -1:
                    continue
                self.ax.axvline(tA)
                
                saveDate = pltRange[0].replace(
                                microsecond=0).isoformat().replace(':', '')
                saveName = '{}_microburst_validation.png'.format(saveDate)
                
                plt.savefig(os.path.join(self.saveDir, saveName))
                self.ax.clear()
 
        return
        
        
    def plot(self, timeA, dosA, timeB, dosB, tRange):
        """
        This method plots the narrow microburst time range from
        both spacecraft, and annotate with L, MLT, Lat, Lon.
        """
        validIdA = np.where((dosA != -1E31) & 
                            (timeA > tRange[0]) & 
                            (timeA < tRange[1]))[0]
        validIdB = np.where((dosB != -1E31) & 
                            (timeB > tRange[0]) & 
                            (timeB < tRange[1]))[0]
        if len(validIdA) == 0 or len(validIdB) == 0:
            return -1

        # Plot the timeseries.
        self.ax.plot(timeA[validIdA], dosA[validIdA], 'r', label='AC-6 A')
        self.ax.plot(timeB[validIdB], dosB[validIdB], 'b', label='AC-6 B')
        #self.bx.plot(timeA[validIdA], self.dataA['Alpha'][validIdA], '--r')
        #self.bx.plot(timeB[validIdB], self.dataB['Alpha'][validIdB], '--b')
        self.ax.set(ylabel='Dos1 [counts/s]', xlabel='UTC',
            title='AC6-{} microburst validation | {}'.format(
                                self.sc_id, tRange[0].date()))
        #self.bx.set(xlabel='UTC', ylabel='Alpha [Deg] (dashed)')
        self.ax.legend(loc=1)

        meanFlag = (np.mean(self.dataA['flag'][validIdA]) + 
            np.mean(self.dataB['flag'][validIdB]))/2
        textStr = ('L={} MLT={}\nlat={} lon={}\ndist={} LCT={}\nflag={}'.format(
            round(np.mean(self.dataA['Lm_OPQ'][validIdA])),
            round(np.mean(self.dataA['MLT_OPQ'][validIdA])),
            round(np.mean(self.dataA['lat'][validIdA])),
            round(np.mean(self.dataA['lon'][validIdA])),
            round(np.mean(self.dataA['Dist_In_Track'][validIdA])),
            round(np.mean(self.dataA['Loss_Cone_Type'][validIdA])),
            round(meanFlag)))
        self.ax.text(0.05, 0.95, textStr, transform=self.ax.transAxes, 
            va='top')
        return 1
        
    def _load_catalog_(self, fPath):
        """
        This method reads in a catalog csv file, and saves it to a dictionary.
        """
        with open(fPath) as f:
            r = csv.reader(f)
            #next(r) # Skip creation header
            keys = next(r)
            rawData = np.array(list(r))
        data = {key:rawData[:, i] for i, key in enumerate(keys)}

        # Now convert the times array(s) to datetime, and all others to float
        for key in itertools.filterfalse(lambda x:'dateTime' in x, 
                                        data.keys()):
                data[key] = data[key].astype(float)
        for key in filter(lambda x: 'dateTime' in x, data.keys()):
            data[key] = np.array([dateutil.parser.parse(i) for 
                                i in data[key]])
        return data 
        
    def _findDateBounds(self):
        """ 
        This function will calculate the unique dates in the catalogue file
        to speed up the plotting reutine.
        """
        dates = np.array([t.date() for t in self.cData['dateTime']])
        uniqueDates = np.array(sorted(set(dates)))
        #print(uniqueDates)

        # Make an array with columns with the unique date, startInd, endInd.
        self.dateArr = np.nan*np.ones((len(uniqueDates), 3), dtype=object)
        # Fill in the self.dateArr with the start/end indicies that correspond
        # to the catalogue.
        for (i, d) in enumerate(uniqueDates):
            #print(dates, d)
            self.dateArr[i, 0] = d
            dateInd = np.where(dates == d)[0]
            #print(d, dateInd)
            self.dateArr[i, 1] = dateInd[0]
            self.dateArr[i, 2] = dateInd[-1]
        return


if __name__ == '__main__':
    sc_id = 'A'
    catPath = ('/home/mike/research/ac6-microburst-scale-sizes/'
                'data/microburst_catalogues/AC6{}_microbursts.txt'.format(
                sc_id))
    
    p = ValidateDetections(sc_id, catPath)
    p.plotLoop()
