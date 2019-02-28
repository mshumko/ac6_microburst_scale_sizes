import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import date
from ac6_microburst_scale_sizes.validation.plot_microbursts import PlotMicrobursts
from matplotlib.widgets import Button

# freqs = np.arange(2, 20, 3)

# fig, ax = plt.subplots()
# plt.subplots_adjust(bottom=0.2)
# t = np.arange(0.0, 1.0, 0.001)
# s = np.sin(2*np.pi*freqs[0]*t)
# l, = plt.plot(t, s, lw=2)

class Browser(PlotMicrobursts):
    def __init__(self, catalog_version, plot_width=5, catalog_save_name=None, width_tol=0.1):
        """
        This class plots the AC6 microbursts and allows the user to browse
        detections in the future and past with buttons. Also there is a button
        to mark the event as a microburst.
        """
        PlotMicrobursts.__init__(self, catalog_version, plot_width=plot_width, 
                                plot_width_flag=False)
        # Filter out events with widths whithin a width_tol.
        if width_tol is not None:
            self.catalog = self.catalog[np.isclose(
                            self.catalog['peak_width_A'], 
                            self.catalog['peak_width_B'], rtol=width_tol)]

        if catalog_save_name is None:
            catalog_save_name = ('AC6_coincident_microbursts_filtered_'
                                'v{}.txt'.format(catalog_version))
        self.current_date = date.min
        self._init_plot()
        #self.lines = [self.ax[0].plot(), self.ax[1].plot()]
        self.index = 0 # Start at row 0 in the dataframe.
        self.filtered_catalog = pd.DataFrame(columns=self.catalog.columns)
        #self.filtered_catalog = pd.DataFrame(data=None, columns=self.catalog.columns, index=self.catalog.index)
        print(self.filtered_catalog)
        self.plot()
        return

    def next(self, event):
        """ Plots the next detection """
        # Just return if at the end of the dataframe.
        if self.index + 1 >= self.catalog.shape[0]:
            return
        self.index += 1
        self.plot()
        return

    def prev(self, event):
        """ Plots the previous detection """
        # Just return if at the end of the dataframe.
        if self.index == 0:
            return
        self.index -= 1
        self.plot()
        return

    def append_microburst(self, event):
        """ 
        Appends the current catalog row to self.filtered_catalog which will then
        be saved to a file for later processing.
        """
        #print(pd.DataFrame(self.catalog.iloc[self.index]))
        if not hasattr(self, 'filtered_catalog'):
            self.filtered_catalog = self.catalog.iloc[self.index].copy()
        print(self.filter_catalog)
        # self.filtered_catalog = self.filtered_catalog.append(self.catalog.iloc[self.index])
        # print(self.filter_catalog)
        return

    def plot(self):
        """ 
        Given a self.current_row in the dataframe, make a space-time plot 
        """
        current_row = self.catalog.iloc[self.index]
        self._clear_ax()
        if current_row['dateTime'].date() != self.current_date:
            print('Loading data from {}...'.format(current_row['dateTime'].date()), 
                    end=' ', flush=True)
            # Load current day AC-6 data if not loaded already
            self.load_ten_hz_data(current_row.dateTime.date())
            self.current_date = current_row.dateTime.date()
            print('done.')
           
        self.make_plot(current_row, savefig=False)
        self.ax[0].set_title('AC6 Microburst Browser\n {} {}'.format(
                        current_row['dateTime'].date(), 
                        current_row['dateTime'].time()))
        self.ax[0].set_ylabel('mean-subtracted dos1rate\n[counts/s]')
        self.ax[1].set_ylabel('mean-subtracted dos1rate\n[counts/s]')
        self.ax[1].set_xlabel('UTC')
        
        self._print_aux_info(current_row)
        plt.draw()
        return

    def _print_aux_info(self, current_row):
        """ Print separation info as well as peak width info to the canvas. """
        self.textbox.clear()
        self.textbox.axis('off')
        s = ('Lag_In_Track = {} s\nDist_In_Track = {} km\n'
                    'Dist_total = {} km\npeak_width_A = {} s\n'
                    'peak_width_B = {} s'.format(
                    round(current_row['Lag_In_Track'], 1), 
                    round(current_row['Dist_In_Track'], 1), 
                    round(current_row['Dist_Total'], 1), 
                    round(current_row['peak_width_A'], 2), 
                    round(current_row['peak_width_B'], 2)))
        self.textbox.text(0, 1, s, va='top')
        return

    def _clear_ax(self):
        [a.clear() for a in self.ax]
        return 

    def _init_plot(self):
        """
        Initialize subplot objects and text box.
        """
        _, self.ax = plt.subplots(2, figsize=(8, 7))
        plt.subplots_adjust(bottom=0.2)
        self.textbox = plt.axes([0.1, 0.05, 0.1, 0.075])
        self.textbox.axis('off')
        return

    def _save_filtered_catalog(self):
        # Remove duplicates

        return


callback = Browser(6)
callback.filter_catalog(filterDict={'Dist_Total':[100, 200]})
axprev = plt.axes([0.59, 0.05, 0.1, 0.075])
axburst = plt.axes([0.7, 0.05, 0.1, 0.075])
axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next)
bprev = Button(axprev, 'Previous')
bprev.on_clicked(callback.prev)
bmicroburst = Button(axburst, 'Microburst')
bmicroburst.on_clicked(callback.append_microburst)
plt.show()