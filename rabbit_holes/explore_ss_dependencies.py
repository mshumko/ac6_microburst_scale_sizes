import numpy as np
import sys
import csv
import matplotlib.pyplot as plt
import dateutil.parser
import os


class ExploreDependencies:
    def __init__(self, cPath):
        self.cat = self._load_catalog(cPath)
        return

    def filter(self, type='flash'):
        """
        This method filters the catalog to only include specific
        burstTypes as catagorized by eye.
        """
        idb = np.where(self.cat['burstType'] == type)[0]

        for key in self.cat.keys():
            self.cat[key] = self.cat[key][idb]
        return

    def _load_catalog(self, cPath):
        """
        This method loads in the coincident microburst catalog.
        """
        with open(cPath, 'r') as f:
            r = csv.reader(f)
            keys = next(r)
            data = np.array(list(r))

        dataDict = {key:data[:, i] for i, key in enumerate(keys)}
        for tkey in filter(lambda x: 'time' in x.lower(), keys):
            dataDict[tkey] = np.array([dateutil.parser.parse(t) for t in dataDict[tkey]])
        for tkey in filter(lambda x: ('time' not in x.lower()) and ('bursttype' not in x.lower()), keys):
            dataDict[tkey] = np.array(list(map(float, dataDict[tkey])))
        return dataDict


if __name__ == '__main__':
    c = ExploreDependencies(os.path.join('/home/mike/research/'
        'ac6-microburst-scale-sizes/data/coincident_microbursts_catalogues', 
        'flash_catalogue_v2_sorted.txt'))
    
