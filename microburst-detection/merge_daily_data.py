# This script will merge files from the daily_microburst_catalogues folder to
# to merged_microburst_catalogues that will have all of the microburst 
# detections from the entire mission. 

import os
import glob
import csv

def mergeDailyFiles(sc_id, inDir, outPath):
    """
    This function goes through text files in inDir that contain the characters
    'AC6<sc_id>', and copies their contents into a file in outPath, excluding 
    the headers.
    """
    # Get all files in directory, sort so files are in temporal order.
    fileNames = sorted(glob.glob(
        os.path.join(inDir, 'AC6{}*.txt'.format(sc_id.upper())) ))
    
    with open(outPath, 'w', newline='') as sF: # Open empty merged file
        writer = csv.writer(sF)
        
        for (i, iFname) in enumerate(fileNames): # Loop over files
            with open(iFname, 'r') as iF:
                reader = csv.reader(iF)
                next(reader)
                if i == 0:
                    writer.writerow(next(reader))
                else:
                    next(reader)
                for line in reader:
                    # Don't copy the header, except for first file.
                    # if ('#' in line[0]) and (i > 0): 
                    #     continue 
                    writer.writerow(line)
    return
    
if __name__ == '__main__':
    for sc_id in ['a', 'b']:
        inPath = ('/home/mike/research/ac6-microburst-scale-sizes/'
            'data/z_daily_microburst_catalogues')
        outPath = ('/home/mike/research/ac6-microburst-scale-sizes/'
            'data/microburst_catalogues/AC6{}_microbursts_v1.txt'.format(sc_id.upper()))
        mergeDailyFiles(sc_id, inPath, outPath)
