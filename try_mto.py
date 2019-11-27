
""" Get image of a galaxy from PanSTARRS1 and run MTO.py"""

import retrieve

from mtolib import utils

import argparse

import numpy as np
from astropy.io import ascii
from astropy.table import Table

import os
import subprocess
from datetime import datetime

#------------------------------------------------

parser = argparse.ArgumentParser(description=
         ''' Downloading images of sky areas with edge-on galaxies \
             from PanSTARRS1 and determining parameters for them.
             \n
             Both methods work: setting coordinates in the command line 
             or getting them from a file. If the coordinates are specified 
             in the command line, the file is ignored. 
             The file should contain a table with the following columns 
             (in the same order): RA, DEC, size.
             \n
             It is assumed that the passband in which the images are loaded,
             and crop factor are the same for every object.
             \n
             Also special parameters for MTO program 
             are the same for every object.
         ''')

parser.add_argument('--ra', help='Right ascension of the centers of galaxies \
                                  (single value or list)',
                    metavar='RA', nargs='*', type=float, default=None)

parser.add_argument('--dec', help='Declination of the centers of galaxies \
                                   (single value or list)',
                    metavar='DEC', nargs='*', type=float, default=None)

parser.add_argument('--size', help='Size of the image to download (in pixels) \
                                    (single value or list)',
                    metavar='SIZE', nargs='*', type=float, default=None)

parser.add_argument('--datafile',help='Name of file containing [RA, DEC, size]',
                    metavar='FILENAME', nargs='?', default=None)

parser.add_argument("--factor", help='Rebinning factor (should be >=1)',
                    metavar='FLOAT', nargs='?', type=float, default=1.) 

parser.add_argument('--image',
                    help='Name of downloaded image fits-file (single value). \
                          In the case of multiple images, all files will be \
                          named "image(No.).fits".',
                    metavar='FILENAME', nargs='?', default='image.fits')

parser.add_argument('--outFile',help='Name of file to write the results to',
                    metavar='FILENAME', nargs='?', default='output.dat')

parser.add_argument("--band", help='Band (filter) name for which galaxy images \
                                    will be downloaded (one of g, r, i, z, y)',
                    metavar='grizy', nargs='?', default='i')

parser.add_argument("--from", dest='start_number', 
                    help='A row number in the Ra-Dec-size table to start from. \
                          Use this only if the previous execution failed, \
                          and you want to continue from a certain row',
                    metavar='NUM', nargs='?', type=int, default=1)


#----------------
# The parameters required for the MTO program will now be described.
#----------------

parser.add_argument('--segImage',
                    help='Name of segmentation image fits-file (single value). \
                          In the case of multiple images, all files will be \
                          named "segImage(No.).fits".',
                    metavar='FILENAME', nargs='?', default='segImage.fits')

parser.add_argument('--csvName',
                    help='Name of resulting MTO csv-file (single value). \
                          In the case of multiple galaxies, all files will be \
                          named "out_mto(No.).csv".',
                    metavar='FILENAME', nargs='?', default='out_mto.csv')

parser.add_argument('--soft_bias', help='Constant to subtract from an image',
                    type=float, default=None)


parser.add_argument('--gain', help='Gain in electrons per ADU', 
                    type=float, default=-1)

parser.add_argument('--bg_mean', help='Mean background',
                    type=float, default=None)

parser.add_argument('--bg_variance', help='Background variance',
                    type=float, default=-1)

parser.add_argument('--alpha', help='Significance level',
                    type=utils.validate_decimal, default=1e-6)

parser.add_argument('--move_factor', help='Moves up the object marker', 
                    type=utils.validate_positive, default=0.5)

parser.add_argument('--min_distance',
                    help='Minimum brightness distance between objects',
                    type=utils.validate_positive, default=0.0)

parser.add_argument('--verbosity', help='Verbosity level (0-2)',
                    type=int, choices=range(0, 3), default=0)


args = parser.parse_args()


#------------------------------------------------
''' In this block, values are simply assigned to the variables. 
    In particular [ra, dec, size] -- depending on how the data
    is received: from the command line or from the file
'''
if None in [args.ra, args.dec, args.size]:
    if args.datafile is None:
        raise argparse.ArgumentTypeError('All parameters [RA,DEC,size] must be \
                                         specified if no data file is provided')
    else:
        dataTable = ascii.read(args.datafile)
        cols = dataTable.colnames # list of the names of columns
        #---------------
        # It is assumed that the table contains the following columns 
        # (in the same order): ra, dec, size
        #---------------
        ra = dataTable[cols[0]]
        dec = dataTable[cols[1]]
        size = dataTable[cols[2]].astype(int)
        number_of_given_objects = len(dataTable)
else:
    # It is assumed that there are many investigated galaxies, 
    # so [ra, dec, size] are lists. 
    # If the coordinates of only one galaxy were specified in the
    # command line, then they will still be converted to lists.
    #---------------
    def convert_to_list(x):
        return x if isinstance(x, (np.ndarray, list, tuple)) else [x]
    #
    # ra = convert_to_list(args.ra)
    # dec = convert_to_list(args.dec)
    # size = convert_to_list(int(float(args.size)))
    ra = args.ra
    dec = args.dec
    size = args.size
    #---------------
    # Checking equality of the length of lists ra, dec, size
    list_sizes = tuple([len(ra), len(dec), len(size)])
    if len(set(list_sizes)) > 1:
        raise argparse.ArgumentTypeError('Lists of RA, DEC, size \
                                          have different sizes')
    number_of_given_objects = len(ra)


factor = args.factor
band = args.band
outputFileName = args.outFile

# In the case of multiple images, all downloaded 
# images will be named "image(No).fits" and saved in the folder "images",
# segmentation images will be named "segImage(No).fits"
# and saved in the folder "segImages"
# CSV-files after running mto.py will have following names: "out_mto(No.).csv"
# and they will be saved in the folder "outCSVs"
#---------------
if number_of_given_objects > 1:
    Ndig = len(str(number_of_given_objects)) # number of digits
    numbers = [ str(n) for n in range(1,number_of_given_objects+1) ]
    #
    if not os.path.exists('images'):
        os.mkdir('images')
    imageNames = ['./images/image{}.fits'.format(n.rjust(Ndig,'0'))
                                                            for n in numbers]
    #
    if not os.path.exists('segImages'):
        os.mkdir('segImages')
    segNames = ['./segImages/segImage{}.fits'.format(n.rjust(Ndig,'0'))
                                                            for n in numbers]
    #
    if not os.path.exists('outCSVs'):
        os.mkdir('outCSVs')
    csvNames = ['./outCSVs/out_mto{}.csv'.format(n.rjust(Ndig,'0'))
                                                            for n in numbers]
else:
    # These are lists (of one item) for uniformity
    imageNames = [args.image] 
    segNames = [args.segImage]
    csvNames = [args.csvName]


#------------------------------------------------

# Getting first image from PanSTARRS1 for testing mto.py
# Running mto.py as test and using output {csv_parameters_demo}
# just to get column names and data format and then delete demo-data

retrieve.main(ra[0], dec[0], output_image=imageNames[0], size=size[0],
                  factor=factor, band=band, format="fits")

subprocess.call('python3 mto.py {} -par_out testCSV.csv'.format(imageNames[0]),
                shell=True)

csv_parameters_demo = Table.read('testCSV.csv')
os.remove('testCSV.csv')


# Creating a resulting table and fill it only by names of columns firstly
out_table = csv_parameters_demo.copy()

# Remove all rows with demo-data after getting names of columns
out_table.remove_rows([i for i in range(len(out_table))]) 

Ncols = len(out_table.colnames)

def convert_Table_row_to_string(table_row):
    ''' The function is needed only to get an astropy.Table row 
        as a string (line) which can be written to a file
    '''
    row_list = list(map(str,table_row.as_void()))
    return ' '.join(row_list)

#------------------------------------------------

# In the case of errors in previous execution, start again from {i_start} row
# If it is the first execution, args.start_number == 0 by default
i_start = args.start_number

if i_start == 1:
    #-----------
    # If {i_start} was > 1, then there would be no need to open the file
    # and write a header, just add new rows to the existing one
    #-----------
    with open(outputFileName, 'w') as out_file:
        out_file.write('# my_id ' + ' '.join(out_table.colnames) + '\n')


# Getting images from PanSTARRS and saving them: module retrieve
for i in range(i_start-1, number_of_given_objects):
    print('---------------------------')
    print()
    print(i+1, datetime.now().strftime('%H:%M:%S'))
    print()
    retrieve.main(ra[i], dec[i], output_image=imageNames[i], size=size[i],
                  factor=factor, band=band, format="fits")
    #-----------
    # Run mto.py via subprocess.call()
    #-----------
    mto_command = []
    mto_command.append('python3 mto.py')
    mto_command.append(imageNames[i])
    mto_command.append('-out ' + segNames[i])
    mto_command.append('-par_out ' + csvNames[i])
    mto_command.append('-gain ' + str(args.gain))
    mto_command.append('-bg_variance ' + str(args.bg_variance))
    mto_command.append('-alpha ' + str(args.alpha))
    mto_command.append('-move_factor ' + str(args.move_factor))
    mto_command.append('-min_distance ' + str(args.min_distance))
    mto_command.append('-verbosity ' + str(args.verbosity))
    if args.soft_bias is not None:
        mto_command.append('-soft_bias ' + str(args.soft_bias))
    if args.bg_mean is not None:
        mto_command.append('-bg_mean ' + str(args.bg_mean))
    #
    mto_run_error_code = subprocess.call(' '.join(mto_command), shell=True)
    #-----------
    if mto_run_error_code != 0:
        # It is possible, so...
        out_table.add_row([-9999 for i in range(len(out_table.colnames))])
        #
        with open(outputFileName, 'a') as out_file:
            out_file.write(' '.join([str(i+1)] + ['-9999']*Ncols) + '\n')
        #
        continue
    #
    #-----------
    # After running MTO: reading csv file to retrieve the computed parameters
    #-----------
    csv_parameters = Table.read(csvNames[i])
    #
    # Find the central galaxy using coordinates of the 
    # image center (size[i]/2, size[i]/2).
    #-----------
    centered_XY = [(csv_parameters['X'][k] - size[i] / factor / 2, \
                    csv_parameters['Y'][k] - size[i] / factor / 2) \
                   for k in range(len(csv_parameters))]
    center_dist = [np.sqrt(point[0]**2 + point[1]**2) for point in centered_XY]
    index_central = list(np.argsort(center_dist))[0]
    #
    out_table.add_row(list(csv_parameters[index_central]))
    #-----------
    # print(list(csv_parameters[index_central]))
    with open(outputFileName, 'a') as out_file:
        out_file.write(str(i+1) + ' ' + \
                       convert_Table_row_to_string(out_table[i-i_start+1])+'\n')

