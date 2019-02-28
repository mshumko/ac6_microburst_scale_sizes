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
        _, self.ax = plt.subplots(2)
        self.lines = [self.ax[0].plot(), self.ax[1].plot()]
        self.current_row = self.catalog.iloc[0] # Start at row 0 in the dataframe.
        self.plot()
        return

    # def next(self, event):
    #     # self.ind += 1
    #     # i = self.ind % len(freqs)
    #     # ydata = np.sin(2*np.pi*freqs[i]*t)
    #     l.set_ydata(ydata)
    #     plt.draw()

    # def prev(self, event):
    #     self.ind -= 1
    #     i = self.ind % len(freqs)
    #     ydata = np.sin(2*np.pi*freqs[i]*t)
    #     l.set_ydata(ydata)
    #     plt.draw()

    def plot(self):
        """ 
        Given a self.current_row in the dataframe, make a space-time plot 
        """
        if self.current_row['dateTime'].date() != self.current_date:
            # Load current day AC-6 data if not loaded already
            self.load_ten_hz_data(self.current_row.dateTime.date())
        self.make_plot(self.current_row)
        return

callback = Browser(6)
callback.filter_catalog(filterDict={'Dist_Total':[0, 25]})
# axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
# axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
# bnext = Button(axnext, 'Next')
# bnext.on_clicked(callback.next)
# bprev = Button(axprev, 'Previous')
# bprev.on_clicked(callback.prev)

plt.show()


########### OLD CODE #############
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.widgets import Button

# freqs = np.arange(2, 20, 3)

# fig, ax = plt.subplots()
# plt.subplots_adjust(bottom=0.2)
# t = np.arange(0.0, 1.0, 0.001)
# s = np.sin(2*np.pi*freqs[0]*t)
# l, = plt.plot(t, s, lw=2)

# class Index:
#     ind = 0

#     def next(self, event):
#         self.ind += 1
#         i = self.ind % len(freqs)
#         ydata = np.sin(2*np.pi*freqs[i]*t)
#         l.set_ydata(ydata)
#         plt.draw()

#     def prev(self, event):
#         self.ind -= 1
#         i = self.ind % len(freqs)
#         ydata = np.sin(2*np.pi*freqs[i]*t)
#         l.set_ydata(ydata)
#         plt.draw()

# callback = Index()
# axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
# axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
# bnext = Button(axnext, 'Next')
# bnext.on_clicked(callback.next)
# bprev = Button(axprev, 'Previous')
# bprev.on_clicked(callback.prev)

plt.show()