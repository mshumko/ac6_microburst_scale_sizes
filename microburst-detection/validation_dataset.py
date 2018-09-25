# This script creates the "true" microburst dataset to validate microburst
# detectors against.
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
from datetime import datetime, timedelta
import dateutil.parser
import numpy as np
#import copy
import sys
import csv
import os

sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

class MakeDataset:
    """
    This class implements the pyplot GUI for the user to identify microbursts
    by hand, and record the times in a file. If a file is not found, this class
    will make one by default. Otherwise it will load the existing dataset and
    add the scatter points to the time series before the user can interract with
    the plot.
    
    The default pyplot GUI commands are the same. To add a microburst detection,
    hover the mouse within ~3 data points of the microburst peak, and click "m".
    This will add a "*" scatter point at the peak location. If you added a point
    by accident, click "r" to remove that point. 
    
    TO navigate, use the "a", "d", "w", and "x" keys to navigate left, right, 
    up, and down, respectifully. If y-axis is in log scale, there will be an 
    error if you press down a few times because a log value of a negative ylim 
    is not defined.
    """
    def __init__(self, date, sc_id, clickWidth=0.5):
        self.sc_id = sc_id
        self.date = date
        self.clickWidth = clickWidth
        self._loadData() # Read in the dosimeter data
        
        fig, self.ax = plt.subplots(figsize=(10, 5))
        self.bx = self.ax.twinx()
        
        # If exists, load the validation dataset 
        self.scatterPts = np.array([])
        if os.path.isfile('validation_dataset.csv'):
            self.validationTimes = self._load_validation_dataset()
            self._plot_existing_detections()
        else:
            self.validationTimes = np.array([])
        return
        
    def plot(self):
        """
        This method plots the dos1rate, dos2rate, and L time series,
        """
        # Plot dos1rate with sigma bounds
        validIdt = np.where(self.d['dos1rate'] != -1E31)[0]
        self.ax.plot(self.d['dateTime'][validIdt], 
                        self.d['dos1rate'][validIdt], 'r', label='dos1rate')
        self.ax.fill_between(self.d['dateTime'][validIdt], 
            self.d['dos1rate'][validIdt]-np.sqrt(self.d['dos1rate'][validIdt]),
            self.d['dos1rate'][validIdt]+np.sqrt(self.d['dos1rate'][validIdt]),
            color='r', alpha=0.5)
        # Plot dos2rate with sigma bounds
        self.ax.plot(self.d['dateTime'][validIdt], 
                        self.d['dos2rate'][validIdt], 'b', label='dos2rate')
        self.ax.fill_between(self.d['dateTime'][validIdt], 
            self.d['dos2rate'][validIdt]-np.sqrt(self.d['dos2rate'][validIdt]),
            self.d['dos2rate'][validIdt]+np.sqrt(self.d['dos2rate'][validIdt]),
            color='b', alpha=0.5)
        # Subplot settings.    
        self.ax.set(yscale='log', xlabel='UTC', ylabel='Dos (counts/s)',
                   ylim=(1, None), 
                   xlim=(self.d['dateTime'][0], self.d['dateTime'][-1])) 
        self.ax.legend(loc=2)
        # Plot L shell on the right-hand side y-axis.
        validL = np.where(self.d['Lm_OPQ'] != -1E31)[0]
        self.bx.plot(self.d['dateTime'][validL], 
                        self.d['Lm_OPQ'][validL], c='k')
        self.bx.set(ylabel='Lm OPQ', ylim=(4, 10))   
        
        # Magical commands to start pyplot's monitoring of key presses.
        self.cid = self.ax.figure.canvas.mpl_connect('key_press_event', self)  
        #self.cid = self.ax.figure.canvas.mpl_connect('button_press_event', self)  
          
        return
        
    def __call__(self, event):
        """
        This method is for the pyplot GUI to record the keypress, and call a
        function that corresponds to the clicked key.
        """
        if event.key not in ['m', 'e', 'a', 'd', 'w', 'x', 's']:
            return

        clickTime = num2date(event.xdata).replace(tzinfo=None) 
        
        # Add or erace detection scatter point.
        if event.key == 'm': self.addPoint(clickTime)
        elif event.key == 'e': self.removePoint(clickTime) # Erace the last detection.

        # Commands to move window left/right/up/down. I did not use "s" key
        # for down since it is the default key to save figure.
        elif event.key == 'a': self.windowLeft()
        elif event.key == 'd': self.windowRight()
        elif event.key == 'w': self.windowUp()
        elif event.key == 'x': self.windowDown()
        
        # Save dataset. Same key will attempt to save the plot.
        # Just say no and it will save.
        elif event.key == 's': self.saveDataset()
        
        self.ax.figure.canvas.draw() # Update plot. 
        return
        
    def show(self):
        plt.show()
        return
        
    def addPoint(self, clickTime):
        """ Find peak around click time, and add to plot """
        dt = timedelta(seconds=self.clickWidth/2)
        validInd = np.where(
                    (self.d['dateTime'] > clickTime - dt) & 
                    (self.d['dateTime'] < clickTime + dt) 
                           )[0]
        imaxC = np.argmax(self.d['dos1rate'][validInd])
        # Add microburst time to self.validationTimes
        self.validationTimes = np.append(
                        self.validationTimes,
                        self.d['dateTime'][imaxC+validInd[0]])
        self.scatterPts = np.append(self.scatterPts, self.ax.scatter(
                        self.d['dateTime'][imaxC+validInd[0]], 
                        self.d['dos1rate'][imaxC+validInd[0]], 
                        s=100, marker='*', color='k'))
        return
        
    def removePoint(self, clickTime):
        """
        Removes the scatter point and time data closest to time 
        of the click
        """
        # Find argument of rthe detection closest to clickTime.
        dts = [(clickTime - t).total_seconds() for t in self.validationTimes]
        idt = np.argmin(np.abs(dts))

        self.validationTimes = np.delete(self.validationTimes, idt)
        # Remove pyplot object (.remove() is a plt function)
        self.scatterPts[idt].remove() 
        # Remove what remains of the removed object (most likeley None)
        self.scatterPts = np.delete(self.scatterPts, idt)
        return 
        
    def windowUp(self):
        """ Move window up. """
        currLims = self.ax.get_ylim()
        dy = currLims[1] - currLims[0]
        newLims = [l + dy for l in currLims]
        self.ax.set_ylim(*newLims)
        return
        
    def windowDown(self):
        """ Move window down """
        currLims = self.ax.get_ylim()
        dy = currLims[1] - currLims[0]
        # Max function help avoid windowing down into neagtive 
        # log values (error).
        newLims = [max([l - dy, 1]) for l in currLims]
        self.ax.set_ylim(*newLims) 
        return
        
    def windowLeft(self):
        """ Move window left """
        currLims = num2date(self.ax.get_xlim())
        dt = currLims[1] - currLims[0]
        newLims = [l - dt for l in currLims]
        self.ax.set_xlim(*newLims)  
        return
            
    def windowRight(self):
        """ Move window right """
        currLims = num2date(self.ax.get_xlim())
        dt = currLims[1] - currLims[0]
        newLims = [l + dt for l in currLims]
        self.ax.set_xlim(*newLims) 
        return
   
    def _plot_existing_detections(self):
        """ 
        In the set-up phase, this method plots the already found microbursts.
        """
        for t in self.validationTimes:
            idt = np.where(self.d['dateTime'] == t)[0]
            self.scatterPts = np.append(self.scatterPts,
                            self.ax.scatter(self.d['dateTime'][idt], 
                            self.d['dos1rate'][idt], 
                            s=100, marker='*', color='k'))
        return

    def _loadData(self):
        """
        Load AC-6 data.

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
    
    def _load_validation_dataset(self):
        """ 
        Open the current dataset (if exists), and populate the plot
        with scatter points.
        """
        with open('validation_dataset.csv', 'r') as f:
            r = csv.reader(f)
            times = np.squeeze(list(r))
        return np.array([dateutil.parser.parse(t) for t in times])
        
    def saveDataset(self, savePath=None):
        """
        This method saves the microburst dataset into a csv file in savePath.
        If savePath is None, it will save to the current directory.
        """
        if savePath is None: savePath = 'validation_dataset.csv'
        print('Saving validation dataset to', savePath)
        
        # Remove duplicates
        saveTimes = num2date(list(set(date2num(self.validationTimes))))
        #saveTimes = [t.replace(tzinfo=None) for t in times]

        with open(savePath, 'w') as f:
            w = csv.writer(f)
            for t in sorted(saveTimes):
                w.writerow([t.replace(tzinfo=None).isoformat()])
        return
        
if __name__ == '__main__':
    sc_id = 'A'
    date = datetime(2016, 10, 14)
    m = MakeDataset(date, sc_id)
    try:
        m.plot()
        m.show()
    except OverflowError:
        print('pyplot overflowed!!')
        pass
    finally:
        m.saveDataset()