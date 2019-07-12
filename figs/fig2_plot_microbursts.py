# This program plots examples of microbursts given times.

from datetime import datetime, timedelta
import matplotlib.pyplot as plt

from mission_tools.ac6.read_ac_data import read_ac_data_wrapper


if __name__ == '__main__':
    plot_width_s = 5
    times = []