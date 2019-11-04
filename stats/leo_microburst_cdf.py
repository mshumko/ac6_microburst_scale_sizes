import os 
import string
import pandas as pd
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

class Microburst_CDF:
    def __init__(self, catalog_version=None, catalog_path=None, max_sep=105):
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
        self.max_sep = max_sep
        self.data_dir = ('/home/mike/research/ac6_microburst_scale_sizes/data')
        self.norm_dir = os.path.join(self.data_dir, 'norm')

        #print(f'Number of microbursts {self.microburst_catalog.shape[0]}')
        return

    def _load_sample_file_(self, path, sum_N=1, offset=0):
        """
        Loads the samples vs. separation CSV file. The sum_N and offset kwargs 
        rebins the sample file. sum_N comnines N separations bins into 1. The
        offset kwarg offsets each bin. The offsets make sense only if they are
        less than sum_N.
        """
        self.samples = pd.read_csv(path, index_col=0)

        if offset != 0:
            self.samples = self.samples[offset:]

        if sum_N > 1 and (not f'{sum_N}km_bins' in path):
            # Resample by taking the sum over sum_N elements and 
            # shift index by offset.
            self.samples = self.samples.groupby(np.arange(self.samples.shape[0])//sum_N).sum()
            # Reasign new index array.
            self.samples = self.samples.set_index(self.samples.index*sum_N+offset)

        self.bin_width = self.samples.index[1] - self.samples.index[0]
        self.sep_bins = self.samples.loc[0:self.max_sep].index
        #print('Sample normalization bins:', self.sep_bins.values)
        return

    def calc_cdf_pdf_stats(self, df, L_lower, L_upper):
        """
        This method calculates the pdf and cdf and errors from a dataframe
        and L shell filtering.
        """
        filtered_catalog = df[(df.Lm_OPQ > L_lower) & (df.Lm_OPQ < L_upper)]
        # Calculate histogram of all events in Dist_Total separation bins
        n, _ = np.histogram(filtered_catalog['Dist_Total'], bins=self.sep_bins)
        # Calculate the weights to scale each element in n by.
        weights = (self.samples['Seconds'].loc[0:self.max_sep-self.bin_width].max()/
                    self.samples['Seconds'].loc[0:self.max_sep-self.bin_width].values)
        n_weighted = np.multiply(weights, n)
        # Calculate the CDF and PDF
        cdf = np.array([np.nansum(n_weighted[i:])/np.nansum(n_weighted) 
                        for i in range(len(n))])
        pdf = np.convolve([-1, 1], cdf, mode='valid')/self.bin_width
        # Calculate the CDF and PDF uncertanties.
        n_prime_std = np.sqrt([sum(n[i:]*weights[i:]**2) for i in range(len(n))])
        cdf_std = cdf*np.sqrt((n_prime_std/sum(n_weighted))**2 + 
                            (n_prime_std[0]/sum(n_weighted))**2)
        pdf_std = np.sqrt(cdf_std[1:]**2 + cdf_std[:-1]**2)/self.bin_width
        return 100*cdf, pdf, 100*cdf_std, pdf_std, filtered_catalog.shape[0]

    def plot_cdf_pdf_all(self, L_array=[4, 5, 6, 7, 8], plot_full_L_range=True, 
                     plot_L=True, sample_name='ac6_norm_all_cdf.csv'):
        """ 
        Plots the CDF and PDF values for the entire radiation belt range and 
        for each L shell bin defined by the L_array argument. 
        """
        _, self.ax = plt.subplots(3, figsize=(8, 8), sharex=True)
        colors=['r', 'b', 'g', 'm', 'c']

        if (plot_L==False) and (plot_full_L_range==False):
            raise  AttributeError("The specified kwargs will not plot anything.")

        # Plot everything
        if plot_full_L_range: self._plot_full_L_range()

        # Plot the L shell bins
        if plot_L:
            for (Li, Lf, color) in zip(L_array[:-1], L_array[1:], colors):
                self._plot_L_bin(Li, Lf, color)

        self._label_and_adjust_subplots()
        return
        
    def save_data(self, path):
        """ Save the CDF/PDF data and errors to a csv file. """
        C, P, C_std, P_std, N = self.calc_cdf_pdf_stats(self.microburst_catalog, 
                                                    4, 8)
        df = pd.DataFrame(data=np.array([C[:-1], P, C_std[:-1], P_std]).T, 
                          columns=['CDF', 'PDF', 'CDF_std', 'PDF_std'], 
                          index=self.sep_bins[:-2])
        df.to_csv(path)
        return

    def _plot_full_L_range(self, dist_file='microburst_cdf_pdf_v4.csv', 
                            norm_name='ac6_norm_all_cdf.csv'):
        """ Plots the CDF and pdf curves for the full L shell range. """
        # Load the file full of number of samples
        self._load_sample_file_(os.path.join(self.norm_dir, norm_name))

        # Load the file containing the pre-calculated cdf and pdf distribution values.
        dist_data = pd.read_csv(os.path.join(self.data_dir, dist_file), index_col=0)

        # Plot the PDF and CDF curves assuming the errors are symmetric about the mean. 
        # The fill_between higlights the 95% CI (2 standard deviations)
        self.ax[0].fill_between(dist_data.index, dist_data['cdf']-2*dist_data['cdf_err'], 
                                            dist_data['cdf']+2*dist_data['cdf_err'], 
                                            facecolor='k', alpha=0.5)
        self.ax[1].fill_between(dist_data.index, dist_data['pdf']-2*dist_data['pdf_err'], 
                                            dist_data['pdf']+2*dist_data['pdf_err'], 
                                            facecolor='k', alpha=0.5)

        self.ax[0].plot(dist_data.index, dist_data['cdf'], c='k', 
                    label=f'all', lw=4, alpha=1)
        self.ax[1].plot(dist_data.index, dist_data['pdf'], c='k', lw=4)
        self.ax[2].step(self.sep_bins, self.samples.loc[:m.max_sep]/10000, c='k')
        return

    def _plot_L_bin(self, Li, Lf, color):
        """ Plots the CDF and pdf curves for one particular L bin. """
        if Li > Lf: Li, Lf = Lf, Li # Reverse order if passed reversed L shells.

        # Load the L-dependent sample file
        self._load_sample_file_(
                os.path.join(self.norm_dir, f'ac6_norm_{Li}_L_{Lf}_5km_bins.csv')
                )
        # Calculate the cdf and pdf values
        C, P, _, _, N = self.calc_cdf_pdf_stats(self.microburst_catalog, Li, Lf)
        # Make CDF, PDF, and sample plots.
        self.ax[0].plot(self.sep_bins[:-1], C, c=color,
                    label=f'{Li} < L < {Lf}')
        self.ax[1].plot(self.sep_bins[:-2], P, c=color)
        self.ax[2].step(self.sep_bins, self.samples.loc[:m.max_sep]/10000, c=color)
        return

    def _label_and_adjust_subplots(self):
        """ 
        Helper function to make the 1000 various adjustments to the subplots 
        """
        self.ax[0].legend()
        self.ax[0].set_xlim(0, 100)
        self.ax[0].set_ylim(bottom=0)
        self.ax[1].set_ylim(bottom=0, top=0.07)
        self.ax[0].set_title('AC6 separation distribution of > 35 keV microbursts in low Earth orbit')
        self.ax[0].set_ylabel('Percent of Microbursts Larger')
        self.ax[1].set_ylabel('Microburst Size Histogram')
        self.ax[2].set_xlabel('AC6 Separation [km]')
        self.ax[2].set_ylabel(r'Samples Per Bin x $10^4$')
        self.ax[2].set_ylim(bottom=0)
        self.ax[2].set_xticks(np.arange(min(self.sep_bins), max(self.sep_bins)+1, 10))

        # Add subplot labels (A), (B), (C)
        for i, ax_i in enumerate(self.ax):
            ax_i.text(0, 1, f'({string.ascii_lowercase[i]})', 
                    transform=ax_i.transAxes, va='top', fontsize=15)
        plt.tight_layout()
        return


if __name__ == "__main__":
    catalog_version = 6
    catalog_path = (f'/home/mike/research/'
                        'ac6_microburst_scale_sizes'
                        '/data/coincident_microbursts_catalogues/'
                        'AC6_coincident_microbursts_sorted'
                        f'_err_v{catalog_version}.txt')
    m = Microburst_CDF(catalog_version=None, catalog_path=catalog_path)
    #_, ax = plt.subplots(3, figsize=(8, 8), sharex=True)
    #m._plot_full_L_range(ax=ax)
    m.plot_cdf_pdf_all()
    # m.plot_cdf_pdf_all(plot_full_L_range=True, plot_L=True, sample_name='ac6_norm_all_1km_bins.csv')
    # m.save_data('/home/mike/research/ac6_microburst_scale_sizes/'
    #             'data/microburst_cdf_pdf_norm_v3.csv')
    plt.show()
