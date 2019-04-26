import os 
import pandas as pd
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

class Microburst_CDF:
    def __init__(self, catalog_version=None, catalog_path=None):
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
        self.max_sep = 100
        print(f'Number of microbursts {self.microburst_catalog.shape[0]}')
        return

    def _load_sample_file_(self, path):
        """
        Loads the samples vs. separation CSV file 
        """
        self.samples = pd.read_csv(path, index_col=0)
        self.bin_width = self.samples.index[1] - self.samples.index[0]
        self.sep_bins = self.samples.loc[0:self.max_sep].index
        return
            
    def n_i(self, di, df, data=None):
        """ 
        Calculates the number of microbursts in data 
        dataframe between separation di and df.
        """
        if data is None:
            data = self.microburst_catalog
        n = sum((data['Dist_Total'] >= di) & 
                (data['Dist_Total'] < df))
        return  n

    def calc_cdf_pdf_stats(self, df, L_lower, L_upper):
        """
        This method calculates the pdf and cdf and errors from a dataframe
        and L shell filtering.
        """
        filtered_catalog = df[(df.Lm_OPQ > L_lower) & (df.Lm_OPQ < L_upper)]
        # Apply the separation bins to the n_i function to get an array.
        n = np.array([self.n_i(bi, bf, filtered_catalog) for bi, bf in 
                    zip(self.sep_bins[:-1], self.sep_bins[1:])]).flatten()
        # Calculate the weights to scale each element in n by.
        weights = (self.samples['Seconds'].loc[0:self.max_sep-self.bin_width].max()/
                    self.samples['Seconds'].loc[0:self.max_sep-self.bin_width].values)
        n_weighted = np.multiply(weights, n)
        # Calculate the CDF and PDF
        cdf = np.array([sum(n_weighted[i:])/sum(n_weighted) for i in range(len(n))])
        pdf = np.convolve([-1, 1], cdf, mode='valid')/self.bin_width
        # Calculate the CDF and PDF uncertanties.
        n_prime_std = np.sqrt([sum(n[i:]*weights[i:]**2) for i in range(len(n))])
        cdf_std = cdf*np.sqrt((n_prime_std/sum(n_weighted))**2 + 
                            (n_prime_std[0]/sum(n_weighted))**2)
        pdf_std = np.sqrt(cdf_std[1:]**2 + cdf_std[:-1]**2)/self.bin_width
        return cdf, pdf, cdf_std, pdf_std, filtered_catalog.shape[0]

    def calc_cdf_pdf_resample(self, df, L_lower, L_upper, N_resample=int(100)):
        """
        This method calculates the pdf and cdf and errors from a dataframe
        and L shell filtering. The CDF and PDF curves are estimated in the 
        same way as self.calc_cdf_pdf_stats, but here we resample the 
        microbursts using a MC method and recalculate the CDF and get 
        uncertanities as a function of separation.
        """
        filtered_catalog = df[(df.Lm_OPQ > L_lower) & (df.Lm_OPQ < L_upper)]
        # Calculate the weights to scale each element in n by.
        weights = (self.samples['Seconds'].loc[0:self.max_sep-self.bin_width].max()/
                    self.samples['Seconds'].loc[0:self.max_sep-self.bin_width].values)

        cdf_array = np.nan*np.ones((len(self.sep_bins)-1, N_resample))
        pdf_array = np.nan*np.ones((len(self.sep_bins)-2, N_resample))

        for i in range(N_resample):
            # Resample the catalog by "flippiping" filtered_catalog.shape[0] number of "coins"
            # one time to get a random sample of microbursts assuming a random false_rate.
            random_flips = np.random.binomial(1, p=1-filtered_catalog.false_rate)
            keep_indicies = filtered_catalog.index[random_flips != 0]
            # Apply the separation bins to the n_i function to get an array.
            n = np.array([self.n_i(bi, bf, filtered_catalog.loc[keep_indicies]) for bi, bf in 
                        zip(self.sep_bins[:-1], self.sep_bins[1:])]).flatten()
            n_weighted = np.multiply(weights, n)
            # Calculate the CDF and PDF
            cdf_array[:, i] = np.array([sum(n_weighted[i:])/sum(n_weighted) for i in range(len(n))])
            pdf_array[:, i] = np.convolve([-1, 1], cdf_array[:, i], mode='valid')/self.bin_width

        # Calculate the CDF and PDF
        cdf = np.mean(cdf_array, axis=1)
        pdf = np.mean(pdf_array, axis=1)
        # Calculate PDF and CDF standard deviations.
        # cdf_err = np.std(cdf_array, axis=1)
        # pdf_err = np.std(pdf_array, axis=1)
        # Calculate the 95% confidence intervals
        cdf_err = np.percentile(cdf_array, [2.5, 97.5], axis=1)
        pdf_err = np.percentile(pdf_array, [2.5, 97.5], axis=1)
        print(cdf_err[1] - cdf_err[0])
        return cdf, pdf, cdf_err, pdf_err, filtered_catalog.shape[0]

    def plot_cdf_pdf(self, L_array=[4, 5, 6, 7, 8], plot_all=True, err_mode='stats', plot_L=False):
        """ Plots the CDF and PDF values. """
        _, ax = plt.subplots(3, figsize=(8, 8), sharex=True)
        c=['r', 'b', 'g', 'm']
        sample_file_dir = ('/home/mike/research/ac6_microburst'
                            '_scale_sizes/data/norm')
        # Running average parameter. Set to 1 to not average.
        n = 1
        if plot_all:
            # Plot the CDF over all L shells in the belts.
            self._load_sample_file_(
                os.path.join(sample_file_dir, 'ac6_norm_all_cdf.csv')
                )
            if err_mode == 'stats':
                C, P, C_err, P_err, N = self.calc_cdf_pdf_stats(self.microburst_catalog, 
                                                        4, 8)
                P = np.convolve(np.ones(n)/n, P, mode='same')
            else: 
                C, P, C_err, P_err, N = self.calc_cdf_pdf_resample(self.microburst_catalog, 
                                                        4, 8)
                P = np.convolve(np.ones(n)/n, P, mode='same')

            # If the errors are symmetric about the mean (standard deviation). 
            if len(C_err.shape) == 1:
                ax[0].fill_between(self.sep_bins[:-1], C-C_err, C+C_err, facecolor='k', 
                                    alpha=0.5)
                ax[1].fill_between(self.sep_bins[:-2], P-P_err, P+P_err, facecolor='k',
                        label=f'4 < L < 8', lw=3, alpha=0.5)
            else:
                ax[0].fill_between(self.sep_bins[:-1], C_err[0], C_err[1], facecolor='k', 
                                    alpha=0.5)
                ax[1].fill_between(self.sep_bins[:-2], P_err[0], P_err[1], facecolor='k',
                        label=f'4 < L < 8', lw=3, alpha=0.5)

            ax[0].plot(self.sep_bins[:-1], C, c='k', 
                        label=f'4 < L < 8 | N = {N}', lw=3, alpha=0.5)
            ax[1].plot(self.sep_bins[:-2], P, c='k', lw=3)
            ax[2].plot(self.sep_bins, self.samples.loc[:m.max_sep]/10000, c='k')
            
        if plot_L:
            for i, (lower_L, upper_L) in enumerate(zip(L_array[:-1], L_array[1:])):
                self._load_sample_file_(
                    os.path.join(sample_file_dir, f'ac6_norm_{lower_L}_L_{upper_L}_5km_bins.csv')
                    )
                if err_mode == 'stats':
                    C, P, C_std, P_std, N = self.calc_cdf_pdf_stats(self.microburst_catalog, 
                                                            lower_L, upper_L)
                else:
                    C, P, C_std, P_std, N = self.calc_cdf_pdf_resample(self.microburst_catalog, 
                                                            4, 8)
                # Running average on the PDF
                P = np.convolve(np.ones(n)/n, P, mode='same')
                ax[0].errorbar(self.sep_bins[:-1], C, c=c[i],
                            label=f'{lower_L} < L < {upper_L} | N = {N}', capsize=5)
                ax[1].errorbar(self.sep_bins[:-2], P, c=c[i], 
                            label=f'{lower_L} < L < {upper_L}')
                ax[2].plot(self.sep_bins, self.samples.loc[:m.max_sep]/10000, c=c[i])
            
        ax[0].legend()
        ax[0].set_xlim(left=0, right=90)
        ax[0].set_ylim(bottom=0)
        ax[1].set_ylim(bottom=0)
        ax[0].set_ylabel('Fraction of Microbursts Larger')
        ax[1].set_ylabel('Microburst PD')
        ax[2].set_xlabel('Separation [km]')
        ax[2].set_ylabel(r'Samples x $10^4$')
        ax[2].set_ylim(bottom=0)
        ax[2].set_xticks(np.arange(min(self.sep_bins), max(self.sep_bins)+1, 10))
        plt.tight_layout()
        plt.show()
        return
        
    def save_data(self, path):
        """ Save the CDF/PDF data and errors to a csv file. """
        C, P, C_std, P_std, N = self.calc_cdf_pdf(self.microburst_catalog, 
                                                    4, 8)
        df = pd.DataFrame(data=np.array([C[:-1], P, C_std[:-1], P_std]).T, 
                          columns=['CDF', 'PDF', 'CDF_std', 'PDF_std'], 
                          index=self.sep_bins[:-2])
        df.to_csv(path)
        return

if __name__ == "__main__":
    catalog_version = 6
    catalog_path = (f'/home/mike/research/'
                        'ac6_microburst_scale_sizes'
                        '/data/coincident_microbursts_catalogues/'
                        'AC6_coincident_microbursts_sorted'
                        f'_err_v{catalog_version}.txt')
    m = Microburst_CDF(catalog_version=None, catalog_path=catalog_path)
    m.plot_cdf_pdf(err_mode='')
    # m.save_data('/home/mike/research/ac6_microburst_scale_sizes/'
    #             'data/microburst_cdf_pdf_norm_v3.csv')
