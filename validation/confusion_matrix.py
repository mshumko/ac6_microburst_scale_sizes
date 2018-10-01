import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date 
from datetime import datetime, timedelta
import dateutil.parser
import sys

import pandas as pd

sys.path.insert(0, '/home/mike/research/'
                    'ac6-microburst-scale-sizes/'
                    'microburst-detection')
import microburst_detection

class ConfusionMatrix(microburst_detection.FindMicrobursts):
    """
    This class uses a validation dataset with a specified microburst 
    detection code to generate a confusion matrix. To run the
    microburst detector, use the getMicroburstIdx() method that is 
    defined in the base class.
    """
    def __init__(self, sc_id, date, 
                validationDataset='validation_dataset.csv'):
        # Load validation data.
        self.vData = pd.read_csv(validationDataset, squeeze=True, 
                                converters={0:dateutil.parser.parse})
        super().__init__(sc_id, date)
        return
        
    def find_microbursts(self, *args, **kwargs):
        """ Wrapper to find microbursts """
        self.getMicroburstIdx(*args, **kwargs)
        return

    def confusionMatrix(self):
        """
        This method will print out the confusion matrix once the 
        microburst detector has been run on the data.
        """
        # Convert detection times to numbers
        detTimes = self.d['dateTime'][self.peakInd]
        # detNum = date2num(detTimes)
        detNum = pd.to_numeric(pd.Series(data=detTimes))
        # Convert validation dataset to numbers
        vNum = pd.to_numeric(self.vData)
        numTrue, self.TPR = self.true_positive_rate(detNum, vNum)
        numFalse, self.FPR = self.false_positive_rate(detNum, vNum)
        print('Number of true microbursts', len(vNum))
        print('Number of detections made', len(self.peakInd), '\n')
        
        print('Total true positive detections', numTrue)
        print('True positive rate', 100*self.TPR, '%\n')
        print('Total false positive detections', numFalse)
        print('False positive rate', 100*self.FPR, '%')
        return

    def true_positive_rate(self, detNum, vNum):
        """ 
        Find the true positive rate == TP/P = True Positive / vNum
        """
        numTrue = sum(np.isin(vNum, detNum))
        TPR = numTrue/len(vNum)
        return numTrue, TPR

    def false_positive_rate(self, detNum, vNum):
        """ 
        Find the false positive rate == FP/N
        FP is the false positive, or false alarm, or Type 1 error
        N is the number of real negative cases in the data
        """
        numFalse = sum(np.logical_not(
                        np.isin(detNum, vNum)))
        N = len(self._pos_cases()) - len(vNum) # Number of real negative cases
        FPR = numFalse/N
        return numFalse, FPR
        
    def plot_detections(self, ax=None):
        """ 
        This method evaluates how well the confusion matrix works 
        """
        if ax is None:
            fig, self.ax = plt.subplots()
        else:
            self.ax = ax
        validInd = np.where(self.d['dos1rate'] != -1E31)[0]
        
        self.ax.plot(self.d['dateTime'][validInd],
                    self.d['dos1rate'][validInd], 'r', label='dos1rate')
        self.ax.plot(self.d['dateTime'][validInd],
                    self.d['dos2rate'][validInd], 'b', label='dos2rate')
        
        # Plot scatter points of valid detections
        vNum = pd.to_numeric(self.vData)
        datNum = pd.to_numeric(pd.Series(data=self.d['dateTime']))
        idv = np.where(np.isin(datNum, vNum))[0]
        self.ax.scatter(self.d['dateTime'][idv],
                    self.d['dos1rate'][idv], c='r', s=100, 
                    label='valid detections')
                    
        # Plot scatter points of detected detections
        #datNum = pd.to_numeric(pd.Series(data=self.d['dateTime']))
        #idv = np.where(np.isin(datNum, vNum))[0]
        self.ax.scatter(self.d['dateTime'][self.peakInd],
                    self.d['dos1rate'][self.peakInd], 
                    c='b', s=100, marker='x',
                    label='detected detections')
                    
        if ax is None: self.ax.legend()
        
        self.ax.set(xlabel='UTC', ylabel='dos1rate [counts/s]', 
                    yscale='log')
        return

    def _pos_cases(self):
        """
        This method is similar to the _checkMicroburstFlag method in
        microburst_detection.py, except it returns the total number
        of data points in that day that could qualify as a detection.
        """
        flag = np.array(self.d['flag'], dtype=int)
        # Identify flags without the ground transmission.
        gTxInd = np.bitwise_and(flag, 1)
        # Identify flags without the cross-spacecraft transmission.
        scTxInd = np.bitwise_and(flag, 2)
        # Index array of bad flags.
        badBool = np.logical_or(gTxInd, scTxInd) 
        goodBool = np.logical_not(badBool)    
        goodInd = np.where(goodBool)[0]
        return goodInd

if __name__ == '__main__':
    dKwargs = {'method':'wavelet', 'thresh':0.1, 'maxWidth':0.25,
                'SIGNIF_LEVEL':0.50}
    c = ConfusionMatrix('A', datetime(2016, 10, 14))
    c.find_microbursts(**dKwargs)
    c.confusionMatrix()
    c.plot_detections()
    plt.show()
