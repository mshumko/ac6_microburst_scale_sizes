import numpy as np
import csv
import dateutil.parser 
import matplotlib.pyplot as plt

class ReplaceErrorVals:
    def __init__(self, sepPath, loadCatPath, saveCatPath=None):
        """
        This class will replace error values in the AC-6 separation and in-track
        time lags in the flash/curtain catalogs.
        """
        self.sepPath = sepPath
        self.loadCatPath = loadCatPath

        if saveCatPath is None:
            self.saveCatPath = loadCatPath
        else:
            self.saveCatPath = saveCatPath
        return

    def loadSeprationFile(self):
        """
        This method will load in the separation file and
        parse the times.
        """
        self.sepDict = {}
        self.sepKeys, self.sepDict = self._load_data(self.sepPath, self.sepDict)

        # Parse times
        # Convert datetimes
        for key in filter(lambda x: 'Date/Time' in x, self.sepKeys):
            self.sepDict[key] = np.array([dateutil.parser.parse(i) 
                                    for i in self.sepDict[key]])
        
        # Convert lags and separations to floats
        for key in filter(lambda x: 'Date/Time' not in x, self.sepKeys):
            self.sepDict[key] = np.array(self.sepDict[key], dtype=float)
        return

    def loadCatalog(self):
        """
        This method reads in the flash or curtain catalog
        """
        self.catDict = {}
        self.catKeys, self.catDict = self._load_data(self.loadCatPath, self.catDict)
        # Parse times
        # Convert datetimes
        for key in filter(lambda x: 'dateTime' in x, self.catKeys):
            self.catDict[key] = np.array([dateutil.parser.parse(i) 
                                    for i in self.catDict[key]])
        for key in filter(lambda x: 'dateTime' not in x and 'burstType' not in x, self.catKeys):
            self.catDict[key] = np.array(self.catDict[key], dtype=float)
        return

    def replaceErrors(self):
        """
        This function looks for error values in the flash/curtain catalog.
        """
        # Find all errorrous values in the catalog file.
        errInd = np.where(self.catDict['Lag_In_Track'] == -1E31)[0]
        #print(errInd)
        # Loop over each error, and find the in-track separation and time 
        # separation values from the separation file given by Bern Blake.
        if 'dateTimeA' in self.catKeys:
            tKey = 'dateTimeA'
        else:
            tKey = 'dateTime'
        for i in errInd:
            catInd = np.where((self.catDict[tKey][i] > self.sepDict['Date/Time'][:-1]) &
                            (self.catDict[tKey][i] < self.sepDict['Date/Time'][1:]) )[0]

            if len(catInd) > 1:
                print('More than one separation time found! \n {}, time={}'.format(catInd,      self.catDict[tKey][i]))
                continue
            elif len(catInd) == 0:
                print('No separation time found! \n {}, time={}'.format(catInd,      self.catDict[tKey][i]))
                continue

            self.catDict['Lag_In_Track'][i] = self.sepDict['Time Separation [sec]'][catInd[0]]
            self.catDict['Dist_In_Track'][i] = self.sepDict['In-Track Separation [km]'][catInd[0]]
        return

    def _load_data(self, path, d):
        """
        This is the basic method that reads in csv data.
        """
        with open(path, 'r') as f:
            raw_data = csv.reader(f)
            keys = next(raw_data)
            data = np.array(list(raw_data), dtype=object)

        # Save data into a dictionary
        for (i, key) in enumerate(keys):
            d[key] = data[:, i]
        return keys, d

    def save_data(self):
        with open(self.saveCatPath, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(self.catKeys) # write header keys
            writer.writerows(zip(*[self.catDict[key] for key in self.catKeys]))

if __name__ == '__main__':
    # r = ReplaceErrorVals('/home/mike/research/ac6/AC6_Separation.csv', 
    #     ('/home/mike/research/ac6-microburst-scale-sizes/'
    #     'data/flash_catalogues/flash_catalogue_v2.txt'))
    for sc_id in ['A', 'B']:
        v = 3
        dName = 'AC6{}_microbursts_v{}.txt'.format(sc_id.upper(), v)
        r = ReplaceErrorVals('/home/mike/research/ac6/AC6_Separation.csv', 
            ('/home/mike/research/ac6-microburst-scale-sizes/'
            'data/microburst_catalogues/{}'.format(dName)))
        r.loadSeprationFile()
        r.loadCatalog()
        r.replaceErrors()
        r.save_data()
    
    # Plot the lifetime AC-6 separation
    # plt.plot(r.sepDict['Date/Time'], r.sepDict['In-Track Separation [km]'])
    # plt.xlabel('Date')
    # plt.ylabel('Separation [km]')
    # plt.title('AC-6 In-Track Separation')
    # plt.tight_layout()
    # plt.show()