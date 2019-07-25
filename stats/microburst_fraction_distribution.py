# This program calculates and plots the fraction of simulatenous to 
# all microburst detections as a function of separation. Uncertanity 
# is first assumed due to Poisson noise and will be then expanded to
# a systematic uncertanity due to transmitter noise.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

class MicroburstFraction:
    def __init__(self, microburst_catalog_name, c_microburst_catalog_name,
                microburst_catalog_dir=None, c_microburst_catalog_dir=None):
        # Specify the coincident catalog directory
        if c_microburst_catalog_dir is None:
            self.c_microburst_catalog_dir = ('/home/mike/research/'
                            'ac6_microburst_scale_sizes/'
                            'data/coincident_microbursts_catalogues')
        else:
            self.c_microburst_catalog_dir = c_microburst_catalog_dir

        # Specify the entire microburst catalog directory
        if microburst_catalog_dir is None:
            self.microburst_catalog_dir = ('/home/mike/research/'
                            'ac6_microburst_scale_sizes/data/'
                            'microburst_catalogues')
        else:
            self.microburst_catalog_dir = microburst_catalog_dir
        
        # Load the catalogs
        self.microburst_catalog = self.load_catalog(
            os.path.join(self.microburst_catalog_dir, microburst_catalog_name))
        self.c_microburst_catalog = self.load_catalog(
            os.path.join(self.c_microburst_catalog_dir, c_microburst_catalog_name))
        return

    def load_catalog(self, path):
        """ 
        Loads a microburst catalog (or) any CSV file with column
        names in header. Path is the absolute path to the catalog
        file. Returns a pandas DataFrame object with the catalog. 
        """
        return pd.read_csv(path)

if __name__ == '__main__':
    microburst_name = 'AC6B_microbursts_v5.txt'
    coincident_catalog_name = 'AC6_coincident_microbursts_sorted_Brady_v6.txt'

    mf = MicroburstFraction(microburst_name, coincident_catalog_name)