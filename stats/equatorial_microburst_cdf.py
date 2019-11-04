import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import os
import string

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
        # print(f'Number of microbursts {self.microburst_catalog.shape[0]}')
        return

    def calc_cdf_pdf(self, df, L_lower, L_upper, bin_width=100):
        """
        This method calculates the pdf and cdf and errors from a dataframe
        and L shell filtering.
        """
        self.bin_width = bin_width
        self.filtered_catalog = df[(df.Lm_OPQ > L_lower) & (df.Lm_OPQ < L_upper)]
        # Map to the magnetic equator
        self.filtered_catalog.loc[:, 'd_equator'] = np.array([
                    self.map2equator(row.lat, row.lon, row.alt, 
                                    row.dateTime, row.Dist_Total) 
                    for _, row in self.filtered_catalog.iterrows()])
        
        # Load the equatorial normalization file.
        self._load_norm(bin_width)

        # Calculate the CDF
        n, _ = np.histogram(self.filtered_catalog['d_equator'], bins=self.norm.index)
        # Before calculating the weights, sum over the appropriate normalization columns.
        L_vals = np.arange(L_lower, L_upper)
        samples = np.zeros_like(self.norm.index)
        for L_col in L_vals:
            samples = samples + self.norm.loc[:, str(float(L_col))]
        # Calculate the weights to scale each element in n by.
        weights = (samples.loc[0:self.norm.index[-1]-self.bin_width].max()/
                    samples.loc[0:self.norm.index[-1]-self.bin_width].values)
        n_weighted = np.multiply(weights, n)

        total_detections = self.filtered_catalog.shape[0]
        cdf = np.array([sum(n_weighted[i:])/sum(n_weighted) for i in range(len(n))])
        pdf = n_weighted/(self.bin_width*sum(n_weighted))

        # Calculate CDF and PDF errors. Assume CDF errors are just due to 
        # Counting statistics like I did in the LEO case.
        cdf_std, pdf_std = self._calc_errors(n, weights)
        return 100*cdf, pdf, 100*cdf_std, pdf_std, total_detections

    def plot_cdf_pdf(self, L_array=[4, 5, 6, 7, 8], plot_all=True):
        """ Plots the CDF and PDF values. """
        _, self.ax = plt.subplots(3, figsize=(8, 8), sharex=True)
        c=['r', 'b', 'g', 'm']
        # sample_file_dir = ('/home/mike/research/ac6_microburst'
        #                     '_scale_sizes/data/norm') 
        if plot_all:
            # Plot the CDF over all L shells in the belts.
            cdf, pdf, cdf_std, pdf_std, N = self.calc_cdf_pdf(
                                                    self.microburst_catalog, 
                                                    4, 8)
            # Plot just the line
            self.ax[0].errorbar(self.norm.index[:-1], cdf, c='k',
                        label=f'all', lw=2, capsize=2)
            self.ax[1].errorbar(self.norm.index[:-1], pdf, c='k', lw=2)
            # Try using fill_between for errors
            self.ax[0].fill_between(self.norm.index[:-1:1], cdf[::1]-2*cdf_std[::1], 
                        cdf[::1]+2*cdf_std[::1], facecolor='k', alpha=0.5)
            self.ax[1].fill_between(self.norm.index[:-1:1], pdf[::1]-2*pdf_std[::1], 
                        pdf[::1]+2*pdf_std[::1], facecolor='k', alpha=0.5)
            # Plot the number of samples.
            all_samples = np.zeros_like(self.norm.index)
            for L_col in L_array[:-1]:
                all_samples = all_samples + self.norm.loc[:, str(float(L_col))]
            self.ax[2].step(self.norm.index, all_samples/1E5, 'k')           

        for i, (lower_L, upper_L) in enumerate(zip(L_array[:-1], L_array[1:])):
            cdf, pdf, _, _, N = self.calc_cdf_pdf(
                                                    self.microburst_catalog, 
                                                    lower_L, upper_L)
            # Plot just the line
            self.ax[0].errorbar(self.norm.index[:-1], cdf, c=c[i],
                        label=f'{lower_L} < L < {upper_L}', capsize=5)
            self.ax[1].errorbar(self.norm.index[:-1], pdf, c=c[i], 
                        label=f'{lower_L} < L < {upper_L}')
            self.ax[2].step(self.norm.index, self.norm[str(float(lower_L))]/1E5, c[i])
        self._label_and_adjust_subplots()
        return

    def _calc_errors(self, n, w, n_trials=10_000):
        """ 
        Given the number of equatorial detections by n for n bins, calculate the PDF
        and CDF errors assuming Poisson statistics. The n array is normalized by an
        array of weights passed as w.
        """
        sqrt_n = np.sqrt(n) # Poisson error from the actual number of observations
        # w*np.sqrt(n) term scales the error by the normalization.
        pdf_std = w*sqrt_n/(self.bin_width*sum(n*w)) # Now normalize it to an actual PDF.

        # Calculate the standard deviation range of n values and calculate the lower 
        # and upper cdf bounds. The cdf_std will then be half of the difference.
        n_upper = w*(n + sqrt_n)
        n_lower = w*(n - sqrt_n)

        cdf_upper = np.array([sum(n_upper[i:]) for i in range(len(n_upper))])/np.sum(n_upper)
        cdf_lower = np.array([sum(n_lower[i:]) for i in range(len(n_lower))])/np.sum(n_lower)
        cdf_std = (cdf_upper-cdf_lower)/2

        return cdf_std, pdf_std

    def _load_norm(self, bin_width):
        """ 
        Load the equatorial normalization file and rebin if the bin_width 
        is not equal to the index difference.
        """
        norm_dir = '/home/mike/research/ac6_microburst_scale_sizes/data/norm'
        norm_name = 'equatorial_norm.csv'
        norm_path = os.path.join(norm_dir, norm_name)
        self.norm = pd.read_csv(norm_path, index_col=0)
        sep_min = self.norm.index.min()
        sep_max = self.norm.index.max()

        if self.norm.index[1] - self.norm.index[0] != bin_width:
            # Now rebin by the bin sizes.
            self.norm = self.norm.groupby(self.norm.index//bin_width).sum()
            # Replace the consecutive indicies with [0, bin_width, 2*bin_width...] 
            self.norm = self.norm.set_index(
                        np.arange(sep_min, sep_max+1, bin_width))
        return

    def _label_and_adjust_subplots(self):
        """ 
        Helper function to make the 1000 various adjustments to the subplots 
        """
        self.ax[0].legend()
        self.ax[0].set_xlim(left=1, right=self.norm.index[-3])
        self.ax[0].set_ylim(bottom=0)
        self.ax[1].set_ylim(bottom=0)
        self.ax[0].set_title('AC6 equatorial separation distribution of > 35 keV microbursts')
        self.ax[0].set_ylabel('Percent of Microbursts Larger')
        self.ax[1].set_ylabel('Microburst Size Histogram')
        self.ax[-1].set_xlabel('AC6 Equatorial Separation [km]')
        self.ax[-1].set_ylabel(r'Samples Per Bin x $10^5$')
        # Add plot labels (A), (B), (C)
        for i, ax_i in enumerate(self.ax):
            ax_i.text(0, 1, f'({string.ascii_lowercase[i]})', 
                    transform=ax_i.transAxes, va='top', fontsize=15)
        plt.tight_layout()
        return

    def deltaLat(self, d, alt):
        """
        Calculate the half of the change in angle for a spacecraft at
        an altitude alt and separated by a distance d.
        """
        return np.rad2deg(np.arcsin(d/(2*(Re+alt))))

    def map2equator(self, lat, lon, alt, time, d):
        """ Maps to magnetic equator assuming d is soly in latitude. """
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
    eq = Microburst_Equatorial_CDF(6)
    eq.plot_cdf_pdf()
    plt.show()