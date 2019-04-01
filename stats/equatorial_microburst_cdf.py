import os 
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

import pandas as pd

import IRBEM

Re = 6371 # Earth radius, km

class Microburst_Equatorial_CDF:
    def __init__(self, catalog_version, catalog_path=None):
        """
        This class calculates the microburst CDF (and PDF) 
        distrbutions. This class essentially takes the 
        relevant code from the microburst_cdf.ipynb and
        generalizes a few parts to make a clean plots of
        the PDF and CDF.

        Load a filtered catalog since this class does not
        do any filtering.
        """
        if catalog_path is None:
            catalog_path = (f'/home/mike/research/'
                        'ac6_microburst_scale_sizes'
                        '/data/coincident_microbursts_catalogues/'
                        'AC6_coincident_microbursts_sorted_'
                        f'Brady_v{catalog_version}.txt')
        # Load catalog.
        self.microburst_catalog = pd.read_csv(catalog_path)
        self.model = IRBEM.MagFields(kext='OPQ77')
        print(f'Number of microbursts {self.microburst_catalog.shape[0]}')
        return


    def calc_cdf_pdf(self, df, L_lower, L_upper):
        """
        This method calculates the pdf and cdf and errors from a dataframe
        and L shell filtering.
        """
        self.filtered_catalog = df[(df.Lm_OPQ > L_lower) & (df.Lm_OPQ < L_upper)]
        # Map to the magnetic equator
        self.filtered_catalog['d_equator'] = np.array([
                    self.map2equator(row.lat, row.lon, row.alt, 
                                    row.dateTime, row.Dist_Total) 
                    for _, row in self.filtered_catalog.iterrows()])

        # for row in self.filtered_catalog.loc[['lat', 'lon', 'alt', 'dateTime', 'Dist_Total']].iterrows():
        #     print(row)
        # self.filtered_catalog['d_equator'] = np.array([
        #             self.map2equator(*row) 
        #             for row in self.filtered_catalog.loc[:, ('lat', 'lon', 'alt', 'dateTime', 'Dist_Total')]])

        # Calculate the CDF
        self.bin_width = 25
        self.bins = np.arange(0, 2000, self.bin_width)
        total_detections = self.filtered_catalog.shape[0]
        cdf = np.array([len(np.where(self.filtered_catalog.d_equator > d)[0])/
                                total_detections for d in self.bins])
        pdf = (cdf[:-1] - cdf[1:])/self.bin_width
        return cdf, pdf

    def plot_cdf_pdf(self, L_array=[4, 5, 6, 7, 8], plot_all=True):
        """ Plots the CDF and PDF values. """
        _, ax = plt.subplots(2, figsize=(8, 8), sharex=True)
        c=['r', 'b', 'g', 'm']
        sample_file_dir = ('/home/mike/research/ac6_microburst'
                            '_scale_sizes/data/norm')            

        for i, (lower_L, upper_L) in enumerate(zip(L_array[:-1], L_array[1:])):
            cdf, pdf = self.calc_cdf_pdf(self.microburst_catalog, 
                                                    lower_L, upper_L)
            ax[0].errorbar(self.bins, cdf, c=c[i], 
                        label=f'{lower_L} < L < {upper_L}', capsize=5)
            ax[1].errorbar(self.bins[:-1], pdf, c=c[i], 
                        label=f'{lower_L} < L < {upper_L}')

        if plot_all:
            # Plot the CDF over all L shells in the belts.
            cdf, pdf = self.calc_cdf_pdf(self.microburst_catalog, 
                                                    4, 8)
            ax[0].errorbar(self.bins, cdf, c='k', 
                        label=f'4 < L < 8', lw=3, capsize=5)
            ax[1].errorbar(self.bins[:-1], pdf, c='k', 
                        label=f'4 < L < 8', lw=3)
            
        ax[0].legend()
        ax[0].set_xlim(left=-1, right=2000)
        ax[0].set_ylim(bottom=0)
        ax[1].set_ylim(bottom=0)
        ax[0].set_ylabel('Microburst fraction')
        ax[1].set_ylabel('Microburst PD')
        ax[1].set_xlabel('Size [km]')
        #ax[2].set_xticks(np.arange(min(self.sep_bins), max(self.sep_bins)+1, 10))
        plt.tight_layout()
        plt.show()
        return

    def deltaLat(self, d, alt):
        """
        Calculate the half of the change in angle for a spacecraft at
        an altitude alt and separated by a distance d.
        """
        return np.rad2deg(np.arcsin(d/(2*(Re+alt))))

    def map2equator(self, lat, lon, alt, time, d):
        """ Maps to magnetix equator assuming d is soly in latitude. """
        # Define the coordinates of the two spacecraft.
        dLat = self.deltaLat(d, alt)
        X1 = {'x1':alt, 'x2':lat-dLat, 'x3':lon, 'dateTime':time}
        X2 = {'x1':alt, 'x2':lat+dLat, 'x3':lon, 'dateTime':time}
        # Run IRBEM
        X1_equator = self.model.find_magequator(X1, None)['XGEO']
        X2_equator = self.model.find_magequator(X2, None)['XGEO']
        # Calculate the separations.
        self.d_equator = Re*np.linalg.norm(X1_equator-X2_equator)
        return self.d_equator

if __name__ == "__main__":
    m = Microburst_Equatorial_CDF(6)
    m.plot_cdf_pdf()