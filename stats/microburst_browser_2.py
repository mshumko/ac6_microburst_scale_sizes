import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import date
from ac6_microburst_scale_sizes.validation.plot_microbursts import PlotMicrobursts
from matplotlib.widgets import Button

catalog_save_dir = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
                    'coincident_microbursts_catalogues')

class Browser(PlotMicrobursts):
    def __init__(self, catalog_version, plot_width=5, 
                catalog_save_name=None, width_tol=0.1):
        """
        This class plots the AC6 microbursts and allows the user to browse
        detections in the future and past with buttons. Also there is a button
        to mark the event as a microburst.
        """
        PlotMicrobursts.__init__(self, catalog_version, plot_width=plot_width, 
                                plot_width_flag=False, make_plt_dir_flag=False)
        # Filter out events with widths whithin a width_tol.
        if width_tol is not None:
            self.catalog = self.catalog[np.isclose(
                            self.catalog['peak_width_A'], 
                            self.catalog['peak_width_B'], rtol=width_tol)]

        if catalog_save_name is None:
            self.catalog_save_name = ('AC6_coincident_microbursts_sorted_'
                                'v{}.txt'.format(catalog_version))
        else:
            self.catalog_save_name = catalog_save_name

        self.current_date = date.min
        self._init_plot()
        self.index = 0 # Start at row 0 in the dataframe.
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
        if not hasattr(self, 'microburst_idx'):
            self.microburst_idx = np.array([self.index])
        else:
            self.microburst_idx = np.append(self.microburst_idx, self.index)
        print('Micorburst saved at', self.catalog.iloc[self.index].dateTime)
        return

    def plot(self):
        """ 
        Given a self.current_row in the dataframe, make a space-time plot 
        """
        print('Index position = {}/{}'.format(self.index, self.catalog.shape[0]-1))
        current_row = self.catalog.iloc[self.index]
        self._clear_ax()

        if current_row['dateTime'].date() != self.current_date:
            # Load current day AC-6 data if not loaded already
            print('Loading data from {}...'.format(current_row['dateTime'].date()), 
                    end=' ', flush=True)
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

    def save_filtered_catalog(self):
        # Remove duplicates
        self.microburst_idx = np.unique(self.microburst_idx)
        save_path = os.path.join(catalog_save_dir, self.catalog_save_name)
        print('Saving filtered catalog to {}'.format(save_path))
        df = self.catalog.iloc[self.microburst_idx]
        df.to_csv(save_path, index=False)
        return


callback = Browser(6)
callback.filter_catalog(filterDict={'Dist_Total':[100, 200]})

# Define button axes.
axprev = plt.axes([0.59, 0.05, 0.1, 0.075])
axburst = plt.axes([0.7, 0.05, 0.1, 0.075])
axnext = plt.axes([0.81, 0.05, 0.1, 0.075])

# Define buttons and their actions.
bnext = Button(axnext, 'Next', hovercolor='g')
bnext.on_clicked(callback.next)
bprev = Button(axprev, 'Previous', hovercolor='g')
bprev.on_clicked(callback.prev)
bmicroburst = Button(axburst, 'Microburst', hovercolor='g')
bmicroburst.on_clicked(callback.append_microburst)

# Initialize the GUI
plt.show()
# Save the catalog.
callback.save_filtered_catalog()