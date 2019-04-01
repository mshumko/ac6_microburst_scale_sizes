import pandas as pd
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

class Microburst_CDF:
    def __init__(catalog_version, catalog_path=None):
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
                        'Brady_v{catalog_version}.txt')
        # Load catalog.
        microburst_catalog = pd.read_csv(catalog_path)
        print(f'Number of microbursts {microburst_catalog.shape[0]}')
        return
            