# This script creates the "true" microburst dataset to validate microburst
# detectors against.
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
from datetime import datetime, timedelta
import numpy as np
#import copy
import sys
import csv

sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

class MakeDataset:
    """
    This class implements the pyplot GUI for the user to identify microbursts
    by hand, and record the times in a file.
    
    winWidth is the width of the view window. Keys 'a' and 'd' move to the 
    previous, and next windows, respectifully.
    """
    def __init__(self, date, sc_id, winWidth=30):
        self.sc_id = sc_id
        self.date = date
        self.winWidth = winWidth
        self._loadData() # Read in the data
        self.trueData = []
        
        fig, self.ax = plt.subplots(figsize=(10, 5))
        self.bx = self.ax.twinx()
        return
        
    def plot(self):
        """
        This method returns the times, and count rates for the next plot.
        """
        validIdt = np.where(self.d['dos1rate'] != -1E31)[0]
        self.ax.plot(self.d['dateTime'][validIdt], 
                        self.d['dos1rate'][validIdt], label='dos1rate')
        self.ax.fill_between(self.d['dateTime'][validIdt], 
            self.d['dos1rate'][validIdt]-np.sqrt(self.d['dos1rate'][validIdt]),
            self.d['dos1rate'][validIdt]+np.sqrt(self.d['dos1rate'][validIdt]),
            color='r', alpha=0.5)
        self.ax.plot(self.d['dateTime'][validIdt], 
                        self.d['dos2rate'][validIdt], label='dos2rate')
        self.ax.fill_between(self.d['dateTime'][validIdt], 
            self.d['dos2rate'][validIdt]-np.sqrt(self.d['dos2rate'][validIdt]),
            self.d['dos2rate'][validIdt]+np.sqrt(self.d['dos2rate'][validIdt]),
            color='r', alpha=0.5)
            
        self.ax.set(yscale='log', xlabel='UTC', ylabel='Dos (counts/s)', ylim=(1, None)) 
        self.ax.legend(loc=2)
        
        validL = np.where(self.d['Lm_OPQ'] != -1E31)[0]
        self.bx.plot(self.d['dateTime'][validL], 
                        self.d['Lm_OPQ'][validL], c='k')
        self.bx.set(ylabel='Lm OPQ', ylim=(4, 10))   
        
        self.cid = self.ax.figure.canvas.mpl_connect('key_press_event', self)  
        self.cid = self.ax.figure.canvas.mpl_connect('button_press_event', self)  
        self.detections = []      
        return
        
    def __call__(self, event):
        """
        This method is for the pyplot GUI to record the keypress.
        """
        #if event.inaxes!=self.ax.axes: return
        #print('click', event)
        ix, iy = event.xdata, event.ydata
        print('Clicked on point', num2date(ix), iy)
        if event.key == 'm':
            # Add microburst time to self.trueData
            self.trueData.append(num2date(ix).replace(tzinfo=None))
            self.lastDet = self.ax.scatter(ix, iy, s=20, marker='*')
        elif event.key == 'r':
            # Remove the last detection.
            self.trueData = self.trueData[:-1]
            self.lastDet.remove()
        elif event.key == 'a':
            # Slide plot window to the left by winWidth
            currLims = num2date(self.ax.get_xlim())
            dt = currLims[1] - currLims[0]
            newLims = [l - dt for l in currLims]
            self.ax.set_xlim(*newLims)  
        elif event.key == 'd':
            # Slide plot window to the left by winWidth
            currLims = num2date(self.ax.get_xlim())
            dt = currLims[1] - currLims[0]
            newLims = [l + dt for l in currLims]
            self.ax.set_xlim(*newLims)     
        return
        
    def show(self):
        plt.show()
        return
        
        
    def _loadData(self):
        """
        Load the AC-6 data.

        10 Hz data keys:
            'alt', 'lat', 'lon', 'dos1l', 'dos1m', 'dos1rate', 'dos2l', 'dos2m',
            'dos2rate', 'dos3l', 'dos3m', 'dos3rate', 'flag', 'Subcom', 'Lm_OPQ', 
            'Bmag_OPQ', 'MLT_OPQ', 'InvLat_OPQ', 'Loss_Cone_Type', 'K_Z', 'Lstar_Z',
            'hmin_Z', 'Alpha', 'Beta', 'Dist_In_Track', 'Lag_In_Track', 
            'Dist_Cross_Track_Horiz', 'Dist_Cross_Track_Vert', 'Dist_Total'
        """
        self.d = read_ac_data.read_ac_data_wrapper(self.sc_id, self.date,
            dType='10Hz', plot=False)
        return
        
    def saveCatalog(self, savePath=None):
        """
        This method saves the microburst catalog into a csv file in savePath.
        If savePath is None, it will save to the current directory.
        """
        if savePath is None: savePath = 'true_catalog.csv'
        
        with open(savePath, 'a') as f:
            w = csv.writer(f)
            for t in self.trueData:
                w.writerow([t.isoformat()])
        return
        
if __name__ == '__main__':
    sc_id = 'A'
    date = datetime(2016, 10, 14)
    m = MakeDataset(date, sc_id)
    m.plot()
    m.show()
    m.saveCatalog()
        
