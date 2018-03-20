# This script will map the scale sizes to the magnetic equator.

import IRBEM
import csv
import numpy as np

Re = 6371 # km
def deltaLat(d, alt):
    """
    This function calculates the latitude angle from an arc length 
    (effectively straight line) given in d, at an altitude alt.
    d and alt must be in units of km.
    """
    dLat = 180/np.pi*d/(Re+alt)
    return dLat

class EquatorScaleSize:
    def __init__(self, cDir):
        """

        """
        return

if __name__ == '__main__':
    ss = EquatorScaleSize()