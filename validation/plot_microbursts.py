# This program plots the detections that make it through into the 
# CDF plots.
import os
import numpy as np
import pandas as pd
from datetime import datetime, date, timedelta
import matplotlib.pyplot as plt

# Paths that are used everywhere in the class
CATALOG_DIR = ('/home/mike/research/ac6_microburst_scale_sizes/data/'
                'coincident_microbursts_catalogues/')
AC6_DATA_PATH = lambda sc_id: ('/home/mike/research/ac6/ac6{}/'
                                'ascii/level2'.format(sc_id))
PLOT_SAVE_DIR = '/home/mike/Desktop/ac6_microburst_validation'

class PlotMicrobursts:
    def __init__(self, catalog_version, plot_width=5, catalog_name=None,
                plot_save_dir=None, plot_width_flag=True, 
                make_plt_dir_flag=True, load_sorted_catalog=True):
        """
        This class plots the detections from the coincident
        microburst catalog with a set of default filters.
        The user can supply other filter values in the 
        filterDict.
        """
        self.plot_width = timedelta(seconds=plot_width)
        self.catalog_name = catalog_name
        if plot_save_dir is None:
            self.plot_save_dir = PLOT_SAVE_DIR
        else:
            self.plot_save_dir = plot_save_dir

        if make_plt_dir_flag and (not os.path.exists(self.plot_save_dir)):
            os.mkdir(self.plot_save_dir)
            print('Made directory:', self.plot_save_dir)
        self.plot_width_flag = plot_width_flag
        self.load_sorted_catalog = load_sorted_catalog
        self.load_catalog(catalog_version)
        return

    def filter_catalog(self, defaultFilter=True, filterDict={}):
        """
        Apply filters to the catalog file. The default filter
        applies a filter that should always be used. The dafault
        filter filters out the events with a temporal_CC < 0.8,
        detections made outside of 4 < L < 8, filters out detections
        that have a spatial_CC-0.3 > temporal_CC. Lastly, the default
        filters removes detections made above the US and SAA to avoid
        noise from ground transmitters and saturation in the SAA.

        Whenever you want to filter other variables, supply them to
        the filterDict dictionary. filterDict's keys are the
        catalog keys and the values is a 2 element list with upper 
        and lower bounds. If you want a one-sided filter e.g. 
        Dist_Total greater than 50 km, pass in 
        filterDict={'Dist_Total':[50, 9999]}.
        """
        if defaultFilter:
            # High CC filter
            self.catalog = self.catalog[self.catalog['time_cc'] >= 0.8]
            # Rad belt filter
            self.catalog = self.catalog[
                                        (np.abs(self.catalog['Lm_OPQ']) > 4) &
                                        (np.abs(self.catalog['Lm_OPQ']) < 8)
                                        ]
            # Curtain filter
            self.catalog = self.catalog[self.catalog['time_cc'] > 
                                        self.catalog['space_cc']+0.3]
            # USA filter
            self.catalog = self.catalog[
                    ((self.catalog['lon'] > -60) | (self.catalog['lon'] < -140))|
                    ((self.catalog['lat'] > 70) | (self.catalog['lat'] < 15))
                                        ]
            # SAA filter
            # self.catalog = self.catalog[
            #         ((self.catalog['lon'] > 30)  | (self.catalog['lon'] < -116)) |
            #         ((self.catalog['lat'] < -90) | (self.catalog['lat'] > 0))
            #                     ]
            # Significane above baseline filter
            #self.catalog = self.catalog[self.catalog['peak_std'] > 3]

            for key, vals in filterDict.items():
                self.catalog = self.catalog[
                    ((self.catalog[key] > min(vals)) & 
                    (self.catalog[key] < max(vals)))
                                ]
        return

    def load_catalog(self, catalog_version):
        """
        Load the csv microburst catalog. Any key that has the word time in
        it is converted to am array of datetime objects. Time keys are:
        dateTime, time_spatial_A, and time_spatial_B.
        """
        # Load a catalog that has been sorted by a person or not.
        if self.catalog_name is None:
            if self.load_sorted_catalog: 
                catalog_name = f'AC6_coincident_microbursts_sorted_v{catalog_version}.txt'
            else:
                catalog_name = f'AC6_coincident_microbursts_v{catalog_version}.txt'
        else:
            catalog_name = self.catalog_name

        catalog_path = os.path.join(CATALOG_DIR, catalog_name)
        self.catalog = pd.read_csv(catalog_path)
        # Convert the catalog times to datetime objects
        for timeKey in ['dateTime', 'time_spatial_A', 'time_spatial_B']:
            self.catalog[timeKey] = pd.to_datetime(self.catalog[timeKey])
        return

    def loop(self, **kwargs):
        """
        Loops over the detections that made it through the filters and 
        make space-time plots.
        """
        current_date = date.min

        _, self.ax = plt.subplots(2)

        for _, row in self.catalog.iterrows():
            if row.dateTime.date() != current_date:
                # Load current day AC-6 data if not loaded already
                self.load_ten_hz_data(row.dateTime.date())
                current_date = row.dateTime.date()
            self.make_plot(row, **kwargs)
            self.ax[0].clear()
            self.ax[1].clear()
        plt.close()
        return

    def load_ten_hz_data(self, day):
        """
        Load the 10 Hz AC6 data from both spacecraft on date.
        """
        time_keys = ['year', 'month', 'day', 'hour', 'minute', 'second']
        dayStr = '{0:%Y%m%d}'.format(day)
        pathA = os.path.join(AC6_DATA_PATH('a'), 
                'AC6-A_{}_L2_10Hz_V03.csv'.format(dayStr))
        pathB = os.path.join(AC6_DATA_PATH('b'), 
                'AC6-B_{}_L2_10Hz_V03.csv'.format(dayStr))
        self.ac6a_data = pd.read_csv(pathA, na_values='-1e+31')
        self.ac6a_data['dateTime'] = pd.to_datetime(self.ac6a_data[time_keys])
        self.ac6b_data = pd.read_csv(pathB, na_values='-1e+31')
        self.ac6b_data['dateTime'] = pd.to_datetime(self.ac6b_data[time_keys])
        return

    def make_plot(self, row, **kwargs):
        """
        This method takes in a dataframe row from the catalog and makes a 
        space/time plot.
        """
        mean_subtracted = kwargs.get('mean_subtracted', True)
        savefig = kwargs.get('savefig', True)
        log_scale = kwargs.get('log_scale', False)
        plot_dos2_and_dos3 = kwargs.get('plot_dos2_and_dos3', True)
        plot_legend = kwargs.get('plot_legend', True)
        ax = kwargs.get('ax', self.ax)
        time_guide_flag = kwargs.get('time_guide_flag', False)

        df_time_a, df_time_b, df_space_a, df_space_b = self._get_filtered_plot_data(row)
        if mean_subtracted:
            df_time_a.loc[:, 'dos1rate'] -= df_time_a.loc[:, 'dos1rate'].mean()
            df_time_b.loc[:, 'dos1rate'] -= df_time_b.loc[:, 'dos1rate'].mean()
            df_space_a.loc[:, 'dos1rate'] -= df_space_a.loc[:, 'dos1rate'].mean()
            df_space_b.loc[:, 'dos1rate'] -= df_space_b.loc[:, 'dos1rate'].mean()

            if plot_dos2_and_dos3:
                df_time_a.loc[:, 'dos2rate'] -= df_time_a.loc[:, 'dos2rate'].mean()
                df_time_b.loc[:, 'dos2rate'] -= df_time_b.loc[:, 'dos2rate'].mean()
                #df_space_a.loc[:, 'dos2rate'] -= df_space_a.loc[:, 'dos2rate'].mean()
                df_time_a.loc[:, 'dos3rate'] -= df_time_a.loc[:, 'dos3rate'].mean()
                #df_space_a.loc[:, 'dos3rate'] -= df_space_a.loc[:, 'dos3rate'].mean()
                
        ax[0].plot(df_time_a['dateTime'], df_time_a['dos1rate'], 'r', label='AC6-A dos1')
        if plot_dos2_and_dos3:
            ax[0].plot(df_time_a['dateTime'], df_time_a['dos2rate'], 'r:', label='AC6-A dos2')
            ax[0].plot(df_time_a['dateTime'], df_time_a['dos3rate'], 'r--', label='AC6-A dos3')
            
        ax[0].plot(df_time_b['dateTime'], df_time_b['dos1rate'], 'b', label='AC6-B')
        if plot_dos2_and_dos3:
            ax[0].plot(df_time_b['dateTime'], df_time_b['dos2rate'], 'b:', label='AC6-B dos2')
        if time_guide_flag:
            ax[0].axvline(row.at['dateTime'])
        if plot_legend:
            ax[0].legend(loc=1)
        ax[1].plot(df_space_a['dateTime'], df_space_a['dos1rate'], 'r', label='AC6-A')
        ax[1].plot(df_space_b['dateTime'], df_space_b['dos1rate'], 'b', label='AC6-B')

        if log_scale:
            ax[0].set_yscale('log')
            ax[1].set_yscale('log')

        # Print peak width if it exists in the catalog.
        if set(['peak_width_A', 'peak_width_B']).issubset(row.index) and self.plot_width_flag:
            s = 'peak_width_A = {} s\npeak_width_B = {} s'.format(
                    round(row['peak_width_A'], 2), round(row['peak_width_B'], 2))
            ax[0].text(0, 1, s, transform=ax[0].transAxes, va='top')
        if savefig:
            save_name = '{0:%Y%m%d_%H%M%S}_ac6_validation_dist_total_{1}.png'.format(
                        row['dateTime'], round(row['Dist_Total']))
            plt.savefig(os.path.join(self.plot_save_dir, save_name))
        return

    def _get_filtered_plot_data(self, row):
        df_time_a = self.ac6a_data[
                            (self.ac6a_data['dateTime'] > row['dateTime']-self.plot_width/2) & 
                            (self.ac6a_data['dateTime'] < row['dateTime']+self.plot_width/2)
                            ]
        df_time_b = self.ac6b_data[
                            (self.ac6b_data['dateTime'] > row['dateTime']-self.plot_width/2) & 
                            (self.ac6b_data['dateTime'] < row['dateTime']+self.plot_width/2)
                            ]
        if row.at['dateTime'] == row.at['time_spatial_A']:
            df_space_a = df_time_a
            df_space_b = self.ac6b_data[
                            (self.ac6b_data['dateTime'] > row['time_spatial_B']-self.plot_width/2) & 
                            (self.ac6b_data['dateTime'] < row['time_spatial_B']+self.plot_width/2)
                            ]
            df_space_b.loc[:, 'dateTime'] -= timedelta(seconds=row.at['Lag_In_Track'])
        elif row.at['dateTime'] == row.at['time_spatial_B']:
            df_space_a = self.ac6a_data[
                            (self.ac6a_data['dateTime'] > row['time_spatial_A']-self.plot_width/2) & 
                            (self.ac6a_data['dateTime'] < row['time_spatial_A']+self.plot_width/2)
                            ]
            df_space_a.loc[:, 'dateTime'] += timedelta(seconds=row.at['Lag_In_Track'])
            df_space_b = df_time_b
        else:
            raise(ValueError('No space matches found!'))
        return df_time_a, df_time_b, df_space_a, df_space_b

if __name__ == '__main__':
    p = PlotMicrobursts(6, catalog_name='AC6_coincident_microbursts_sorted_Brady_v6.txt')
    p.filter_catalog(filterDict={'Dist_Total':[60, 70]})
    # p.catalog = p.catalog[np.isclose(p.catalog['peak_width_A'], 
    #                       p.catalog['peak_width_B'], rtol=0.1)]
    p.loop(mean_subtracted=True)

    # p = PlotMicrobursts(6)
    # p.filter_catalog(filterDict={'Dist_Total':[60, 200]})
    # p.catalog = p.catalog[np.isclose(p.catalog['peak_width_A'], 
    #                       p.catalog['peak_width_B'], rtol=0.1)]
    # p.loop()