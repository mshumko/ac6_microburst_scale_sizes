# This program plots examples of microbursts given times.

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates
import numpy as np
import string

from mission_tools.ac6.read_ac_data import read_ac_data_wrapper
from ac6_microburst_scale_sizes.validation.plot_microbursts import PlotMicrobursts

plt.rcParams.update({'font.size': 15})

class PlotExamples(PlotMicrobursts):
    def __init__(self, catalog_version, plot_width, t0_times):
        """
        This class is a child of the PlotMicrobursts class and makes 
        time-aligned and space-aligned plots.
        """
        PlotMicrobursts.__init__(self, catalog_version, plot_width=plot_width, 
                                plot_width_flag=False, make_plt_dir_flag=False)
        self.t0_times = t0_times
        return

    def plot_examples(self):
        """ 
        This method
        """
        # Get the rows from the catalog that have the variables necessary 
        # for plotting.
        rows = [None]*len(self.t0_times)
        for i in range(len(rows)):
            idx = np.where(self.t0_times[i] == self.catalog.dateTime)[0]
            # Check that one unique time was found.
            if len(idx) != 1:
                raise ValueError('None or multiple indicies found! Check t0_times array')
            rows[i] = self.catalog.iloc[idx[0]]
        # Initialize plotting environment
        self._plot_setup()

        for i, row in enumerate(rows):
            # Loop over t0_times.
            # First load AC6 10Hz data for that day
            self.load_ten_hz_data(row.dateTime.date())
            # Make plots for that day.
            self.make_plot(row, savefig=False, ax=self.ax[:, i], plot_legend=False,
                            mean_subtracted=False, plot_dos2_and_dos3=False)
            self.ax[1, i].xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=2))
            self.ax[1, i].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))

            # Add text to each subplot
            # Separation info
            self.ax[0, i].text(0.99, 0.99, f's = {abs(int(round(row.Dist_In_Track)))} km', 
                            transform=self.ax[0, i].transAxes, va='top', ha='right', fontsize=15)
            self.ax[1, i].text(0.99, 0.99, f'dt = {abs(int(round(row.Lag_In_Track)))} s',
                            transform=self.ax[1, i].transAxes, va='top', ha='right', fontsize=15)
            
        plt.show()
        return

    def _plot_setup(self):
        """
        Helper method to set up the plotting environment.
        """
        self.fig, self.ax = plt.subplots(2, len(self.t0_times), figsize=(12, 6))

        for i in range(len(self.t0_times)):
            self.ax[0, i].get_xaxis().set_visible(False)

        # Set up plot labels.
        self.fig.text(0.5, 0.01, 'UTC', ha='center', va='center')
        self.fig.text(0.015, 0.5, 'dos1 [counts/s]', ha='center', 
                    va='center', rotation='vertical')

        # subplot titles
        for i in range(len(self.t0_times)):
            self.ax[0, i].set_title(f'{self.t0_times[i].date()}')

        plt.subplots_adjust(left=0.07, right=0.99, hspace=0.1)

        # subplot labels
        for i in range(len(self.t0_times)):
            self.ax[0, i].text(0, 0.99, f'({string.ascii_letters[2*i]})', va='top',
                                transform=self.ax[0, i].transAxes, fontsize=20)
            self.ax[1, i].text(0, 0.99, f'({string.ascii_letters[2*i+1]})', va='top',
                                transform=self.ax[1, i].transAxes, fontsize=20)
        return

if __name__ == '__main__':
    plot_width_s = 5
    t0_times = [
                datetime(2016, 9, 30, 1, 56, 29, 800000),
                datetime(2016, 11, 23, 4, 46, 46, 600000),
                datetime(2017, 4, 23, 13, 44, 43, 600000)
                ]
    p = PlotExamples(6, plot_width_s, t0_times)
    p.plot_examples()