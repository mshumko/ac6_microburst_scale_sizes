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
    def __init__(self, catalog_version, plot_width=5, catalog_save_name=None):
        PlotMicrobursts.__init__(self, catalog_version, plot_width=plot_width)

        if catalog_save_name is None:
            catalog_save_name = ('AC6_coincident_microbursts_filtered_'
                                'v{}.txt'.format(catalog_version))
        self.current_date = date.min
        # Subplot objects
        _, self.ax = plt.subplots(2, figsize=(8, 7))
        plt.subplots_adjust(bottom=0.2)
        self.lines = [self.ax[0].plot(), self.ax[1].plot()]
        self.index = 0 # Start at row 0 in the dataframe.
        #self.current_row = self.catalog.iloc[self.index] 
        self.plot()
        return

    def next(self, event):
        """ Plots the next detection """
        # Just return if at the end of the dataframe.
        if self.index >= self.catalog.shape[0]:
            return
        self.index += 1
        #self.current_row = self.catalog.iloc[self.index]
        #self._clear_ax()
        self.plot()
        plt.draw()
        return

    def prev(self, event):
        """ Plots the previous detection """
        # Just return if at the end of the dataframe.
        if self.index == 0:
            return
        self.index -= 1
        #self.current_row = self.catalog.iloc[self.index]
        #self._clear_ax()
        self.plot()
        plt.draw()

    def plot(self):
        """ 
        Given a self.current_row in the dataframe, make a space-time plot 
        """
        current_row = self.catalog.iloc[self.index]
        if current_row['dateTime'].date() != self.current_date:
            # Load current day AC-6 data if not loaded already
            self.load_ten_hz_data(current_row.dateTime.date())
            self.current_date = current_row.dateTime.date()
        self._clear_ax() # Clear anything left over from prior plot.
        self.make_plot(current_row)
        self.ax[0].set_title('AC6 Microburst Browser\n {} {}'.format(
                        current_row['dateTime'].date(), 
                        current_row['dateTime'].time()))
        self.ax[0].set_ylabel('mean-subtracted dos1rate\n[counts/s]')
        self.ax[1].set_ylabel('mean-subtracted dos1rate\n[counts/s]')
        self.ax[1].set_xlabel('UTC')
        # self.ax[1].text()
        return

    def _clear_ax(self):
        [a.clear() for a in self.ax]
        return 

callback = Browser(6)
callback.filter_catalog(filterDict={'Dist_Total':[0, 25]})
axprev = plt.axes([0.59, 0.05, 0.1, 0.075])
axburst = plt.axes([0.7, 0.05, 0.1, 0.075])
axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next)
bprev = Button(axprev, 'Previous')
bprev.on_clicked(callback.prev)
bmicroburst = Button(axburst, 'Microburst')
plt.show()