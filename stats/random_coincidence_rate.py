# This program calculates the random microburst rate for coincident 
# microbursts.

import numpy as np
import os
# import matplotlib.pyplot as plt

import pandas as pd

class Random_Coincidence:
    def __init__(self, coincident_name, microburst_a_name, microburst_b_name, output_path):
        """
        Given a dataset from input_path, calculate the random 
        microburst coicidence rate for each coincident microburst
        identified in input_path.
        """
        self.coincident_microbursts = self._load_catalog(coincident_name)
        return

    def _load_catalog(self, path):
        """ Loads the microburst catalog """
        catalog = pd.read_csv(path)
        catalog['dateTime'] = pd.to_datetime(catalog['dateTime'])
        return catalog