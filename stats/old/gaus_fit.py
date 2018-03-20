import numpy as np
import scipy.optimize

class GausFit:
    def __init__(self, BASELINE_SUBTRACT=False):
        """
        This code will fit a timeseries with a single Gaussian
        for now, and may be extended to multiple Gaussians 
        later on.
        """
        self.BASELINE_SUBTRACT = BASELINE_SUBTRACT
        return

    def fitData(self, x, t, p0=None):
        """
        This function will fit a given timeseries by a gaussian
        """
        if self.BASELINE_SUBTRACT:
            raise NotImplementedError('Baseline subtraction not implemented yet')

        if p0 is None:
            # Generic microburst at the center with 50 ms sigma and max A
            p0 = [np.max(x), t[len(t)//2], 50E-3] 
        # Fit the data
        popt, pcov = scipy.optimize.curve_fit(self.nGaus, t, x, p0=p0, 
            sigma=np.sqrt(x)) # Sigma will be wrong if it was baseline subtracted.
        perr = np.sqrt(np.diag(pcov))
        return popt, perr

    def nGaus(self, t, *args):
        """
        NAME:    nGaus(self, t, *args)
        USE:     N-peaked Gaussian function
        INPUTS:  independent variable t, and *args are the Gaussian parameters
                 with the format p11, p12, p13, p21, p22, p23... where the 
                 first indicie is the peak number, the second indicie 
                 represents the (amplitude, peak center, and sigma) for each Gaussian.
        RETURNS: Numpy array of amplitude values.
        AUTHOR:  Mykhaylo Shumko
        MOD:     2018-02-14
        """
        # Unpack arg tuple if args are passed as an array.
        if len(args) == 1: args = args[0]    
        assert len(args) % 3 == 0, ('Invalid number of arguments. '
                                    'Need 3 arguments for each Gaussian.')
        val = 0.0
        args = np.array(args)
        
        #  This for loop sums over all of the gaussians.
        for i in range(0,len(args),3):
            x = np.divide(np.power((t-args[i+1]), 2),
                (2*np.power(args[i+2], 2)))
            val += args[i]*np.exp(-x.astype(float))
        return val