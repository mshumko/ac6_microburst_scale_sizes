# This is GUI code for visializing detections made from an AC-6 flash 
# database.
 
import tkinter
import tkinter.filedialog
from tkinter import ttk
import os
import sys
import numpy as np
import csv
import itertools
import dateutil.parser
from datetime import datetime, timedelta
import time

# Plotting libraries
import matplotlib.lines as mlines
import matplotlib.dates
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.gridspec as gridspec

# Personal libraries to plot
sys.path.append('/home/mike/research/mission-tools/ac6/')
import read_ac_data
 
class MicroburstBrowser(ttk.Frame):
    """The adders gui and functions."""
    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.root = parent
        self.style = ttk.Style()
        self.style.theme_use('clam')
        # Default time width to plot
        self.tWidth = tkinter.IntVar(self.root)
        self.tWidth.set(5)
        # Current plot date (to avoid constantly reloading data.)
        self.currentDate = datetime.min

        self.init_gui()
        return

    def init_gui(self):
        """Builds GUI."""
        self.root.title('Micoroburst Data Browser')
        self.root.geometry('1075x775')
        self.root.option_add('*tearOff', 'FALSE') # for menu items.
 
        self.grid(column=6, row=10, sticky='nsew') # Overall grid for widgets.
        #self.grid(column=5, row=10, sticky=tkinter.N+tkinter.S+tkinter.E+tkinter.W)
 
        # Set the menubar
        self.menubar = tkinter.Menu(self.root)
 
        self.menu_file = tkinter.Menu(self.menubar)
        self.menu_file.add_command(label='Open', command=self.load_database)
        self.menu_file.add_command(label='Save', command=self.save_database)
        self.menu_file.add_command(label='Exit', command=self._quit)
 
        self.menu_edit = tkinter.Menu(self.menubar)
 
        self.menubar.add_cascade(menu=self.menu_file, label='File')
        #self.menubar.add_cascade(menu=self.menu_edit, label='Edit')
 
        self.root.config(menu=self.menubar)

        # The following line is necessary to return the terminal back to user when 
        # the user presses the 'x' button to kill the program.
        self.root.wm_protocol("WM_DELETE_WINDOW", lambda: self._quit())
 
        # Labels and plot canvas that remain constant throughout execution.
        self.database_str = tkinter.StringVar()
        self.database_str.set('Open database:')
        self.database_label = ttk.Label(self.root, 
            textvariable=self.database_str, font=("Helvetica", 12))
        self.database_label.grid(column=0, row=0, columnspan=2, sticky='w')

        self._init_microburst_buttons() # Buttons to catagorize the detecttions.

        # Initialize plots and the microburst list (empty until user selects 
        # a catalog to view).
        self._init_plots()
        self._init_microburst_browser() 

        # Enable capturing keyboard events
        self.root.bind("<Key>", self._keyPress)
        self.root.bind("<Up>", self._changeMicroburst)
        self.root.bind("<Down>", self._changeMicroburst)
        return
 
    def load_database(self):
        """
        This function will open a file browser to load a csv database file.
        """
        # Home dir: os.path.expanduser('~')
        try:
            self.root.database_name =  tkinter.filedialog.askopenfilename(
                initialdir='/home/mike/research/ac6-microburst-scale-sizes/data/flash_catalogues', 
                title="Select database",
                filetypes=(
                    ("txt files","*.txt"),
                    ("csv files","*.csv"), 
                    ("all files","*.*"))
                    )
        except AttributeError as err:
            if "'tuple' object has no attribute 'rfind'" in str(err):
                self.root.database_name = 'none'
                pass
            else:
                raise            
        # ttk.Label(self.root, text='Open database: {}'.format(os.path.basename(
        #     self.root.database_name))).grid(column=0, row=0, columnspan=2)
        self.database_str.set('Open database: {}'.format(os.path.basename(
             self.root.database_name)))
        if self.root.database_name is not 'none':
            self._open_database()

            # Populate the listBox.
            for (t, bt) in zip(self.dataDict['dateTimeA'], self.dataDict['burstType']):
                self.list.insert(tkinter.END, t)
                if bt != 'nan':
                    self.list.itemconfig(tkinter.END, {'bg':'yellow'})
        return

    def plot_timeseries(self, ax):
        """
        This method will handle plotting the origional and shifted AC-6
        time series from a specific energy channel (default to "dos1")
        """
        self.pltTime = dateutil.parser.parse(self.microburstTime)
        if self.currentDate.date() != self.pltTime.date():
            try:
                print('Loading data... hang on.')
                self.dataA = read_ac_data.read_ac_data_wrapper('A', self.pltTime)
                self.dataB = read_ac_data.read_ac_data_wrapper('B', self.pltTime)
                self.currentDate = self.pltTime
            except AssertionError as err:
                print(str(err))
                return

        # Plot the unajusted timeseries
        # Filter out invalid data
        validIdA = np.where((self.dataA['dos1rate'] != -1E31))[0]
        validIdB = np.where((self.dataB['dos1rate'] != -1E31))[0]
        ax.plot(self.dataA['dateTime'][validIdA], 
            self.dataA['dos1rate'][validIdA], 'r', label='AC-6 A')
        ax.plot(self.dataB['dateTime'][validIdB], 
            self.dataB['dos1rate'][validIdB], 'b', label='AC-6 B')      
        tRange = (self.pltTime-timedelta(seconds=self.tWidth.get()/2), 
                  self.pltTime+timedelta(seconds=self.tWidth.get()/2))
        ax.set_xlim(tRange)

        # Auto set ylims
        self._recalcYlim(ax, tRange)
        # Plot the shifted timeseries
        #self.plot_shifted_times(self.ax[0][1])

        # Format times to show just the minutes and seconds
        tfmt = matplotlib.dates.DateFormatter('%M:%S')
        ax.xaxis.set_major_formatter(tfmt)
        return

    def plot_shifted_times(self, ax):
        """
        This method will handle plotting the shifted AC-6
        time series from a specific energy channel (default to "dos1")
        """
        # Look up time shift and separation
        self.idt = np.where(self.dataDict['dateTimeA'] == self.pltTime)[0]
        assert len(self.idt) == 1, 'No matched microburst found in catalog.'
        lag_in_track = self.dataDict['Lag_In_Track'][self.idt[0]]
        dist_in_track = self.dataDict['Dist_In_Track'][self.idt[0]]
        # Shift data
        shiftedTimes = np.array([t + timedelta(seconds=lag_in_track) 
                                for t in self.dataA['dateTime']])
        # Filter out bad data
        validIdA = np.where((self.dataA['dos1rate'] != -1E31))[0]
        validIdB = np.where((self.dataB['dos1rate'] != -1E31))[0]
        
        ax.clear()
        ax.text(0.05, 0.95, 'In track lag={} s ({} km)'.format(
                round(lag_in_track, 1), dist_in_track),
                fontsize=12, va='top', ha='left', transform=ax.transAxes)
        # Add AC6 legend
        ac6a = mlines.Line2D([], [], color='red',
                          markersize=15, label='AC-6 A')
        ac6b = mlines.Line2D([], [], color='blue',
                          markersize=15, label='AC-6 B')
        ax.legend(handles=[ac6a, ac6b], loc=1)
        # Plot time series
        ax.plot(shiftedTimes[validIdA], 
            self.dataA['dos1rate'][validIdA], 'r', label='AC-6 A')
        ax.plot(self.dataB['dateTime'][validIdB], 
            self.dataB['dos1rate'][validIdB], 'b', label='AC-6 B') 
        # Adjust ax limits and format time.
        tRange = (self.pltTime-timedelta(seconds=self.tWidth.get()/2), 
                  self.pltTime+timedelta(seconds=self.tWidth.get()/2))
        ax.set_xlim(tRange)
        ax.set(xlabel='UTC [MM:SS]', ylabel='dos1 [counts/s]')
        self._recalcYlim(ax, tRange)
        tfmt = matplotlib.dates.DateFormatter('%M:%S')
        ax.xaxis.set_major_formatter(tfmt)
        return

    def save_database(self):
        """
        This method will save the catalog with the sorted detections.
        """
        f = tkinter.filedialog.asksaveasfile(mode='w', defaultextension=".txt", 
            initialdir='/home/mike/research/ac6-microburst-scale-sizes/data/flash_catalogues',
            initialfile=os.path.basename(self.root.database_name).replace('.', '_sorted.'), 
            filetypes=(("txt files","*.txt"),
                        ("csv files","*.csv"), 
                        ("all files","*.*")) 
            )
        if f:
            keys = np.append(self.keys, 'burstType')
            writer = csv.writer(f)
            writer.writerow(keys) # write header keys
            writer.writerows(zip(*[self.dataDict[key] for key in keys]))
        return
        
    def _recalcYlim(self, ax, tRange):
        """
        This method recalculates the ylimits of the data being 
        plotted since a set_xlim has been applied to ax.

        It is necessary since matplotlib autoscales the plot 
        according to the entire data set loaded, not just what is
        being plotted.
        """
        validIdtA = np.where((self.dataA['dateTime'] > tRange[0]) & 
                            (self.dataA['dateTime'] < tRange[1]) & 
                            (self.dataA['dos1rate'] !=-1E31))[0]
        validIdtB = np.where((self.dataB['dateTime'] > tRange[0]) & 
                            (self.dataB['dateTime'] < tRange[1]) & 
                            (self.dataB['dos1rate'] !=-1E31))[0]
        ymin = np.min(np.concatenate((self.dataA['dos1rate'][validIdtA], 
                            self.dataB['dos1rate'][validIdtB])))
        ymax = np.max(np.concatenate((self.dataA['dos1rate'][validIdtA], 
                            self.dataB['dos1rate'][validIdtB])))
        ax.set_ylim((0.8*ymin, 1.2*ymax))
        return 

    def _quit(self):
        """Exits program."""
        self.root.quit()
        self.root.destroy()
        return

    def _open_database(self):
        """
        This function will load in the microburst catalogues
        from both spacecraft
        """
        data = self._inputDataSkelletons(self.root.database_name)

        with open(self.root.database_name) as f:
            reader = csv.reader(f)
            #next(reader) # Skip header
            next(reader)

            for i, line in enumerate(reader): # Read in the data
                data[i] = line
            # Now format data array into a dictionary for each spacecraft
            self.dataDict = {}
            for i, key in enumerate(self.keys):
                self.dataDict[key] = data[:, i]
            # Loop over all keys but datetimes and set the array types to be floats.
            for key in itertools.filterfalse(lambda x: 'dateTime' in x or 'burstType' in x, self.keys):
                self.dataDict[key] = self.dataDict[key].astype(float)

            # Convert datetimes
            for key in filter(lambda x: 'dateTime' in x, self.keys):
                self.dataDict[key] = np.array([dateutil.parser.parse(i) for i in self.dataDict[key]])
        
        # Create a sort array
        if 'burstType' not in self.dataDict.keys():
            self.dataDict['burstType'] = np.nan*np.ones(len(self.dataDict['dateTimeA']), dtype=object)
            print('created burstType array in dataDict')

        return

    def _inputDataSkelletons(self, fPath):
        """
        Create empty dictionary and arrays for microburst catalogues. 
        """
        # Scan the file to get number of lines.
        with open(fPath, 'r') as f:
            N = sum(1 for row in f) - 1

        # Get data keys and data from file
        with open(fPath, 'r') as f:
            reader = csv.reader(f)
            #next(reader) # Skip first line of header
            self.keys = next(reader)
            # Remove comment and empty space chars
            #self.keys[0] = self.keys[0][2:] 
            # Empty data file
            data = np.nan*np.ones((N, len(self.keys)), dtype=object)
        return data

    def _init_plots(self):
        """
        This method initializes the plotting area (figure and canvas to plot
        the AC-6 data on.)
        """
        f, self.ax = plt.subplots(2, 1, figsize=(8, 7))
        
        # Add canvas
        self.canvas = FigureCanvasTkAgg(f, master=self.root)
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=9, columnspan=5, sticky='nwse')
        self.canvas.show()
        # Add subplot labels
        self.ax[0].set(xlabel='UTC [MM:SS]', ylabel='dos1 [counts/s]', title='None')
        self.ax[1].set(xlabel='UTC [MM:SS]', ylabel='dos1 [counts/s]')
        plt.tight_layout()
        # Add plot controls 
        navToolbar = tkinter.Frame(self.root)
        toolbar = NavigationToolbar2TkAgg(self.canvas, navToolbar)
        navToolbar.grid(column=1, row=10, columnspan=2, sticky='w')
        toolbar.update()
        # Add a spinbox for the user to modify the plot width
        timeWidthSpinBox = tkinter.Frame(self.root)
        timeWidthSpinBox.grid(column=4, row=10, columnspan=2, sticky='e')
        # Spinbox label.
        tkinter.Label(timeWidthSpinBox, text='Plot time width [s]', font=("Helvetica", 10)).grid(column=0, row=0)
        # Spinbox
        sb = tkinter.Spinbox(timeWidthSpinBox, from_=1, to=10, width=10,
            textvariable=self.tWidth, font=("Helvetica", 12)) 
        sb.grid(column=1, row=0)
        return

    def _init_microburst_browser(self):
        """
        This method creates an empty list in column 0 that will be
        populated by microburst times when a catalog is loaded. This
        list is also binded by <<ListBoxSelect>> to execute the 
        self._change_microburst_plot() method when a new microburst is
        selected to plot it.
        """
        listFrame = tkinter.Frame(self.root)
        listFrame.rowconfigure(0, weight=1) # Need to strech frame to bottom
        listFrame.grid(row=1, column=0, sticky='nws', rowspan=10)
        # Listbox that is binded to change microbursts
        self.list = tkinter.Listbox(listFrame,
            font=("Helvetica", 12))
        self.list.grid(row=0, column=0, sticky='nws', rowspan=10)
        self.list.bind('<<ListboxSelect>>', self._change_microburst_plot)
        # Scrollbar
        scrollbar = tkinter.Scrollbar(listFrame, orient="vertical")
        scrollbar.grid(row=0, column=1, rowspan=10, sticky='ns')
        self.list.config(yscrollcommand=scrollbar.set)
        scrollbar.config(command=self.list.yview)
        return

    def _change_microburst_plot(self, event):
        """
        This method handles the program when the user changes the 
        microburst from the list on the left side.
        """
        w = event.widget
        # try, except block in case user clicks microburst 
        # list without loading database first.
        try: 
            index = int(w.curselection()[0])
            self.microburstTime = w.get(index)
        except (UnboundLocalError, IndexError):
            return
        print('You selected microburst %d: "%s"' % (index, self.microburstTime))

        # Plot stuff
        self.ax[0].set_title(self.microburstTime)
        self.plot_timeseries(self.ax[0])
        self.plot_shifted_times(self.ax[1])
        self.canvas.show()  

        # Handle the radio button logic and update database here
        if isinstance(self.dataDict['burstType'][self.idt[0]], str):
            self.burstType.set(self.dataDict['burstType'][self.idt[0]])
        else:
            self.burstType.set('None')
        return

    def _changeMicroburst(self, event):
        """
        This method gets called when the up/down key is pressed and the 
        previous/next microburst will be loaded and plotted.
        """
        print('Clicked ', event.keysym)
        return

    def _keyPress(self, event):
        print("pressed", repr(event.char))
        if repr(event.char) == 'f':
            self.dataDict['burstType'][self.idt[0]] = 'flash'
            self.burstType.set('flash')
        elif repr(event.char) == 'c':
            self.dataDict['burstType'][self.idt[0]] = 'curtain'
            self.burstType.set('curtain')
        elif repr(event.char) == 'n':
            self.dataDict['burstType'][self.idt[0]] = 'neither'
            self.burstType.set('neither')
        elif repr(event.char) == 'a':
            self.dataDict['burstType'][self.idt[0]] = 'ambiguous'
            self.burstType.set('ambiguous')
        return

    def _init_microburst_buttons(self):
        """
        This function intializes the radio buttons to select what type of 
        microburst the user observed and binds them to the appropriate 
        functions to modify the loaded catalog.
        """
        self.burstType = tkinter.StringVar()

        radioFrame = tkinter.Frame(self.root)
        radioFrame.grid(column=4, row=0, columnspan=2, sticky='e')

        self.flash_button = ttk.Radiobutton(radioFrame, text='Flash [f]', 
                                            variable=self.burstType, 
                                            val='flash', 
                                            command=self._update_burst_type)
        self.curtain_button = ttk.Radiobutton(radioFrame, text='Curtain [c]', 
                                            variable=self.burstType, 
                                            val='curtain', 
                                            command=self._update_burst_type)
        self.neither_button = ttk.Radiobutton(radioFrame, text='Neither [n]', 
                                            variable=self.burstType, val='neither', 
                                            command=self._update_burst_type)
        self.ambiguous_button = ttk.Radiobutton(radioFrame, text='Ambiguous [a]', 
                                            variable=self.burstType, val='ambiguous',
                                            command=self._update_burst_type)
        self.flash_button.grid(row=0, column=0, sticky='w')
        self.curtain_button.grid(row=0, column=1, sticky='e')
        self.neither_button.grid(row=0, column=2, sticky='e')
        self.ambiguous_button.grid(row=0, column=3, sticky='e')
        return

    def _update_burst_type(self):
        """
        This function will update the burst type array with whatever radiobutton
        the user chose.
        """
        print(self.burstType.get())
        self.dataDict['burstType'][self.idt[0]] = self.burstType.get()
        return
 
if __name__ == '__main__':
    root = tkinter.Tk()
    MicroburstBrowser(root)
    root.mainloop()
