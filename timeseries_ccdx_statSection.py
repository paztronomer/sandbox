''' Code to get basic statistics for a section of a list of exposures 
'''

import os
import time
import numpy as np
import pandas as pd
import fitsio
import pickle
import logging
import argparse
try:
    import partial
except:
    pass
import multiprocessing as mp
import easyaccess as ea

# Setup logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO,
)
# Dictionary taken from pixcorrect/decaminfo.py
ccdnums = {
    'S29':  1,
    'S30':  2,
    'S31':  3,
    'S25':  4,
    'S26':  5,
    'S27':  6,
    'S28':  7,
    'S20':  8,
    'S21':  9,
    'S22':  10,
    'S23':  11,
    'S24':  12,
    'S14':  13,
    'S15':  14,
    'S16':  15,
    'S17':  16,
    'S18':  17,
    'S19':  18,
    'S8':  19,
    'S9':  20,
    'S10':  21,
    'S11':  22,
    'S12':  23,
    'S13':  24,
    'S1':  25,
    'S2':  26,
    'S3':  27,
    'S4':  28,
    'S5':  29,
    'S6':  30,
    'S7':  31,
    'N1':  32,
    'N2':  33,
    'N3':  34,
    'N4':  35,
    'N5':  36,
    'N6':  37,
    'N7':  38,
    'N8':  39,
    'N9':  40,
    'N10':  41,
    'N11':  42,
    'N12':  43,
    'N13':  44,
    'N14':  45,
    'N15':  46,
    'N16':  47,
    'N17':  48,
    'N18':  49,
    'N19':  50,
    'N20':  51,
    'N21':  52,
    'N22':  53,
    'N23':  54,
    'N24':  55,
    'N25':  56,
    'N26':  57,
    'N27':  58,
    'N28':  59,
    'N29':  60,
    'N30':  61,
    'N31':  62
}

def main_aux(pathlist=None, ccd=None, coo=None, raw=None, extens=0):
    if (pathlist is not None):
        # Load the table, with no constraint on filetype, in case path
        # is extremely large
        dt = [('nite', 'i8'), ('expnum', 'i8'), ('band', '|S10'), 
              ('exptime', 'f8'), ('path','|S200'),]
        tab = np.genfromtxt(
            pathlist,
            dtype=dt,
            comments='#',
            delimiter=',',
            missing_values=np.nan,
        )
    # Remove duplicates
    uarr, uidx, uinverse, ucounts= np.unique(
        tab['path'],
        return_index=True,
        return_inverse=True,
        return_counts=True,
    )
    if (uidx.size != tab.size):
        ndup = tab.size - uidx.size
        logging.info('Removing {0} duplicates'.format(ndup))
        tab = tab[uidx]
    # Get CCD name
    def f(dic, x):
        for a, b in dic.iteritems():
            if b == x:
                return a
    ccdname = f(ccdnums, ccd)
    # If raw, extension is going to be the ccd name
    if raw:
        extens = ccdname
    # Construct the lists to open each image
    for b in np.unique(tab['band'][np.where(tab['band'] != -9999)]):
        s = tab[np.where(tab['band'] == b)]
        N = s.size
        z = zip(s['path'], [coo] * N, [extens] * N, [raw] * N, 
                s['nite'], s['expnum'], s['exptime'])
        #
        # Call the function that calculates the stats, adds exposure info (and
        # obviously, reads the section beforehand)
     
    # [fpazch@deslogin sandbox]$ python timeseries_ccdx_statSection.py object_exposures_fullpath_Y5.csv --ccd 46 --coo 725 2565 865 4096 --raw
    exp=None
    exit()
    #
    # NOTE
    # The goal is to read the same section, from a bunch of images, in  
    # parallel. Then stack the set of section and calculate the median 
    # of each pixel.
    #
    # coo list contains [x0, y0, x1, y1]
    coo = [0, 0, 256, 256]
    # coo: coordinates, ext: extension to read from the FITS file, 
    # delta: side size (pixel) of the square used to sample the FITS file 
    kw_section = {
        'coo' : coo,
        'ext' : 0,
        'delta' : 256,
    }
    
    nproc = mp.cpu_count()
    logging.info('Launch {0} parallel processes'.format(nproc))
    P1 = mp.Pool(processes=nproc)
    
    # chunk = int(4096 / kw_section['delta'] * 2048 / kw_section['delta'])
    # The value for chunk will help prevent memory errors. "Chops the iterable 
    # into a number of chunks which it submits to the process pool as separate 
    # tasks"
    chunk = 4

    try:
        print('CHECK FOR LIST_AUX IN FITS')
        partial_aux = partial(fits_section, **kw_section)
        t0 = time.time()
        # map_async does not block the processes, and executes non-ordered
        box_i = P1.map_async(partial_aux, fpath, chunk)
        t1 = time.time()
    except:
        logging.warning('No available module: partial')
        aux_list = [(fnm, 
                     kw_section['coo'], 
                     kw_section['ext'], 
                     kw_section['delta']) for fnm in fpath]
        t0 = time.time()
        boxi = P1.map_async(fits_section, aux_list, chunk)
        t1 = time.time()
    boxi.wait()
    aux_boxi = boxi.get()
    print(aux_boxi[-1])
    # =================

    # Go pixel by pixel doing the median

def fits_section(x_list):
    ''' Function to read section
    '''
    fname, coo, ext, raw = x_list
    x0, y0, x1, y1 = coo
    if os.path.exists(fname):
        fits = fitsio.FITS(fname)
        # sq = np.copy(fits[ext][y0:y1 , x0:x1])
        # fits.close()
        # return sq
    else:
        logging.error('File {0} does not exists'.format(fname))
        return False
    # For raw exposures, fetermine the offset
    if raw:
        # Read the header info
        hdr = fitsio.read_header(fname, ext)
        dsec = hdr['DATASEC'].strip().strip('[').strip(']').replace(':', ',')
        dsec = map(int, dsec.split(','))
        if (dsec[2] > dsec[3]) and (dsec[0] < dsec[1]):
            logging.info('{0} inverse y-reading direction'.format(fname))
            sec = np.copy(fits[ext][y0 + dsec[3] : y1 + 1, 
                                    x0 + dsec[0] : x1 + 1]) 
        elif (dsec[2] > dsec[3]) and (dsec[0] > dsec[1]):
            logging.info('{0} inverse y and x-reading direction'.format(fname))
            sec = np.copy(fits[ext][y0 + dsec[3] : y1 + 1,
                                    x0 + dsec[1] : x1 + 1])
        elif (dsec[2] < dsec[3]) and (dsec[0] > dsec[1]):
            logging.info('{0} inverse x-reading direction'.format(fname))
            sec = np.copy(fits[ext][y0 + dsec[2] : y1 + 1,
                                    x0 + dsec[1] : x1 + 1])
        elif (dsec[2] < dsec[3]) and (dsec[0] < dsec[1]):
            logging.info('{0} increasing reading direction'.format(fname))
            sec = np.copy(fits[ext][y0 + dsec[2] : y1 + 1,
                                    x0 + dsec[0] : x1 + 1])
        fits.close()
        return sec
    else:
        sec = np.copy(fits[ext][y0 : y1 + 1, x0 : x1 + 1])
        fits.close()
        return sec

if __name__ == '__main__':
    
    hgral = 'Time Series constructor. Calculates basic statistics for a'
    hgral += ' section of the CCD, for either raw or processed images'
    ecl = argparse.ArgumentParser(description=hgral)
    h0 = 'Table of night, expnum, band, exptime, and path to images'
    h0 += ' for which stats should be calculated. Please use column names:'
    h0 += ' nite, expnum, band, exptime, path.'
    h0 += 'Format is space-separated or csv table'
    ecl.add_argument('path', help=h0)
    h1 = 'CCD number on which operate'
    ecl.add_argument('--ccd', help=h1, type=int)
    h2 = 'Coordinates of the section on which calculate statistics. Format:'
    h2 += ' x0 y0 x1 y1'
    ecl.add_argument('--coo', help=h2, type=int, nargs=4)
    h3 = 'Flag. Use if inputs are raw image with overscan'
    ecl.add_argument('--raw', help=h2, action='store_true')
    # Parser
    ecl = ecl.parse_args()

    main_aux(pathlist=ecl.path, ccd=ecl.ccd, coo=ecl.coo, raw=ecl.raw)

    print(ecl)

