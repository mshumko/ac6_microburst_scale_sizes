# This program plots examples of microbursts given times.

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates
plt.rcParams.update({'font.size': 16})
import numpy as np

from mission_tools.ac6.read_ac_data import read_ac_data_wrapper
from ac6_microburst_scale_sizes.validation.plot_microbursts import PlotMicrobursts

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
        #plt.tight_layout(pad=0.6)
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
        self.fig.text(0.5, 0.04, 'UTC', ha='center', va='center')
        self.fig.text(0.06, 0.5, 'dos1rate [counts/s]', ha='center', 
                    va='center', rotation='vertical')

        # subplot titles
        for i in range(len(self.t0_times)):
            self.ax[0, i].set_title(f'{self.t0_times[i].date()}')

        plt.subplots_adjust(left=0.2)
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