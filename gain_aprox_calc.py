""" Simple script to get an aproximate value of the GAIN, based on bias 
and dflats

This is a ad-hoc code, so do not expect it to work with other set of files
without modify it first

http://spiff.rit.edu/classes/phys445/lectures/gain/gain.html
http://www.public.asu.edu/~rjansen/ast598/ast598_jansen2014.pdf

2 methods:

1) Janesick's method
G = mean(ADU, box) / [stdev(ADU, box)]^2

If we use both bias and dome flats:
    (<dflat_1> + <dflat_2>) - (<bias_1> + <bias_2>)
G = -----------------------------------------------
    (std(diff dflat)^2 - std(diff bias)^2) 

2) Have dome flats with different exposure times. Seelct boxes on different
images and do a plot: mean (ADU, box) vs variance (ADU, box)

"""

import os
import glob
import time
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
import fitsio


def get_info01(in_list):
    """ Auxiliary method to get information from FITS files
    Inputs
    - bi: dataframe containing paths to bias
    - df: dataframe containing paths to dflats
    - nite: night
    Returns
    - list of lists, each of the nested containing: nite, expnum, ccd, band, 
    gain_a, gain_b, fits 
    """
    bi, df, nite, ext = in_list
    print('MP in nite:{0}'.format(nite))
    # Get all files, avoid first
    res_nite = []
    # Bias
    for idx, row in bi.iloc[1:].iterrows():
        regx = os.path.join(row['path'], '*fits')
        f = glob.glob(regx)
        for fits in f:
            hdu = fitsio.FITS(fits)
            # Header
            hdr = hdu[ext].read_header()
            nite = int(hdr['NITE'])
            expnum = hdr['EXPNUM']
            ccd = hdr['CCDNUM']
            band = np.nan
            gain_a = hdr['GAINA']
            gain_b = hdr['GAINB']
            ftype = 'bias'
            res_nite.append(
                [ftype, nite, expnum, ccd, band, gain_a, gain_b, fits]
            )
            hdu.close()
    # Dome flats
    for idx, xx in df.iloc[1:].iterrows():
        regx = os.path.join(xx['path'], '*fits')
        l = glob.glob(regx)
        for fits in l:
            hdu = fitsio.FITS(fits)
            # Header
            hdr = hdu[ext].read_header()
            nite = int(hdr['NITE'])
            expnum = hdr['EXPNUM']
            ccd = hdr['CCDNUM']
            band = hdr['BAND'].strip()
            gain_a = hdr['GAINA']
            gain_b = hdr['GAINB']
            ftype = 'dflat'
            res_nite.append(
                [ftype, nite, expnum, ccd, band, gain_a, gain_b, fits]
            )
            hdu.close()
    return res_nite

def eq01(inlist):
    """ Equation for gain
    Inputs
    - fnm_f1, fnm_f2: filenames of both dome flats
    - fnm_b1, fnm_b2: filenames of both bias
    """
    fnm_f1, fnm_f2, fnm_b1, fnm_b2, band, ccd, nite = inlist 
    box_side=256 
    ext=0
    w_amp=1024
    h_amp=2048
    #
    print(band, ccd)
    # Open fits
    h_f1 = fitsio.FITS(fnm_f1)
    h_f2 = fitsio.FITS(fnm_f2)
    h_b1 = fitsio.FITS(fnm_b1)
    h_b2 = fitsio.FITS(fnm_b2)
    # Get header information
    # Note as all the images are from the same CCD, only need to read one 
    # header to setup the region for statistics
    hdr = h_f1[ext].read_header()
    ampA = hdr['DATASECA'].replace('[', '').replace(']', '').replace(':', ',')
    ampA = [int(x) for x in ampA.split(',')]
    ampB = hdr['DATASECB'].replace('[', '').replace(']', '').replace(':', ',')
    ampB = [int(x) for x in ampB.split(',')]
    # Get coordinates of the box for the statistics
    mid_xa = ampA[1] - (w_amp / 2)
    mid_ya = ampA[3] - (h_amp / 2)
    mid_xb = ampB[1] - (w_amp / 2)
    mid_yb = ampB[3] - (h_amp / 2)
    y0a, y1a = mid_ya - (box_side / 2), mid_ya + (box_side / 2) 
    x0a, x1a = mid_xa - (box_side / 2), mid_xa + (box_side / 2)
    y0b, y1b = mid_yb - (box_side / 2), mid_yb + (box_side / 2) 
    x0b, x1b = mid_xb - (box_side / 2), mid_xb + (box_side / 2)
    # Get data, do calculation per amplifier
    f1a = h_f1[ext][y0a : y1a, x0a : x1a]
    f2a = h_f2[ext][y0a : y1a, x0a : x1a]
    b1a = h_b1[ext][y0a : y1a, x0a : x1a]
    b2a = h_b2[ext][y0a : y1a, x0a : x1a]
    #
    f1b = h_f1[ext][y0b : y1b, x0b : x1b]
    f2b = h_f2[ext][y0b : y1b, x0b : x1b]
    b1b = h_b1[ext][y0b : y1b, x0b : x1b]
    b2b = h_b2[ext][y0b : y1b, x0b : x1b]
    # Equation
    def g(f1, f2, b1, b2):
        g = (np.mean(f1) + np.mean(f2)) - (np.mean(b1))
        g /= (np.power(np.std(f2 - f1), 2) - np.power(np.std(b2 - b1), 2))
        return g
    # Gain
    ga = g(f1a, f2a, b1a, b2a)
    gb = g(f1b, f2b, b1b, b2b)
    # Information to be returned
    res = [[nite, band, ccd, 'A', ga], [nite, band, ccd, 'B', gb]]
    # Close HDU
    h_f1.close()
    h_f2.close()
    h_b1.close()
    h_b2.close()
    return res

def method01(bias, dflat, 
             fill_info=False, 
             outname='header_info_tab.csv',
             gainout='calc_gain_perAmp.csv',):
    """ Calculates gain changes  
    Inputs
    - tab: dataframe, must contain columns path, reqnum, unitname (=nite)
    """
    # Define pairs of bias, domeflats. Use bias simply by exposure number 
    # order. 
    # Do calculation by band, by CCD, by amplifier
    # Params to read FITS section
    ext = 0
    a1, a2, a3, a4 = 2047 - 128, 2047 + 128, 511 - 128, 511 + 128 
    b1, b2, b3, b4 = 2047 - 128, 2047 + 128, 1535 - 128, 1535 + 128 
    # By reqnum
    nite = bias['unitname'].unique()
    # Lists to harbor
    bias_data = []
    dflat_data = []
    #
    # Get information for each CCD
    if fill_info:
        # Create a list with sections by nite, to be fed to a mp.Pool
        byNite = []
        for n in nite:
            # Avoid the first bias and first domeflat. For this, sort using
            # path then avoid the first
            #
            # Bias
            bi = bias.loc[bias['unitname'] == n]
            bi = bi.sort(['path']) # sort_values for newr pandas versions
            bi = bi.reset_index(drop=True)
            #
            # Dome flats
            df = dflat.loc[dflat['unitname'] == n]
            df = df.sort(['path'])
            df = df.reset_index(drop=True)
            #
            byNite.append([bi, df, n, ext])
        # Setup Pool
        P1 = mp.Pool(processes=mp.cpu_count())
        res = P1.map_async(get_info01, byNite)
        res.wait()
        res.successful()
        res = res.get()
        # Each element of res is a list of lists [[...], [...], [...]] and we
        # need to decrease nesting by 1 level
        res_flatten = [lev2 for lev1 in res for lev2 in lev1]
        # Save results in dataframe
        info_tab = pd.DataFrame(
            res_flatten,
            columns=['ftype', 'nite', 'exp', 'ccd', 'band', 'ga', 'gb', 'fnm']
        )
        info_tab.to_csv(outname, header=True, index=False)
    else:
        info_tab = pd.read_csv(outname)
    # Do the difference in pairs.
    param = []
    for nx in info_tab['nite'].unique():
        # Bands present in that night, for dome flats   
        bx = info_tab.loc[info_tab['nite'] == nx, 'band'].unique()
        # As we have stored float (nan) and strings, use list comprehension
        # instead of np.isnan
        bx = filter(lambda x: isinstance(x, basestring), bx)
        for bnd in bx:
            cx = info_tab.loc[
                (info_tab['nite'] == nx) & (info_tab['band'] == bnd),
                'ccd'
            ]
            cx = cx.unique()
            for ccd in cx:
                # At this level, select the bias and dflats. Make pairs
                # Bias
                c1 = info_tab['nite'] == nx
                c2 = info_tab['band'] == bnd
                c3 = info_tab['ccd'] == ccd
                c4 = info_tab['ftype'] == 'bias'
                c5 = info_tab['ftype'] == 'dflat'
                bias = info_tab.loc[c1 & c3 & c4]
                dflat = info_tab.loc[c1 & c2 & c3 & c5]
                # Sort byt the filename, to have sequential images
                bias = bias.sort('fnm')
                dflat = dflat.sort('fnm')
                # Reset indices to use them to get data
                bias.reset_index(drop=True, inplace=True)
                dflat.reset_index(drop=True, inplace=True)
                # How many pairs of dflats can be done per subset?
                N1 = np.floor(len(bias.index) / 2.)
                N2 = np.floor(len(dflat.index) / 2.)
                Npair = int( np.min([N1, N2]) )
                # Create the pairs
                for i in range(Npair)[::2]:
                    f_bias1 = bias.iloc[i]['fnm']
                    f_bias2 = bias.iloc[i + 1]['fnm']
                    f_dflat1 = dflat.iloc[i]['fnm']
                    f_dflat2 = dflat.iloc[i + 1]['fnm']
                    # Call the calculation
                    param.append([
                        f_dflat1,
                        f_dflat2, 
                        f_bias1,
                        f_bias2,
                        bnd, 
                        ccd,
                        nx,
                    ])
    # Run in parallel
    P1 = mp.Pool(processes=mp.cpu_count())
    g_ls = P1.map_async(eq01, param)
    g_ls.wait()
    g_ls.successful()
    g_ls = g_ls.get()
    # Decrease une level of nesting to save as dataframe
    g_ls = [x for sublist in g_ls for x in sublist]
    # Save into a pandas
    cols = ['nite', 'band', 'ccd', 'amp', 'gain']
    res = pd.DataFrame(data=g_ls, columns=cols)
    # Write out
    res.to_csv(gainout, index=False, header=True)
    print('Wrriten {0}'.format(gainout))
    return True     
    
    """
        # Do the operations on this night, per CCD
        i_bi = pd.DataFrame(bias_data, 
                            columns=['nite', 'exp', 'ccd', 'ga', 'gb', 'fnm'])
        i_f = pd.DataFrame(
            dflat_data,
            columns=['nite', 'exp', 'ccd', 'band', 'ga', 'gb', 'fnm']
        )
        #
        # Do it for 4 pairs of measurements
        # Do it per band
        for bx in i_f['band'].unique():
            # Do it per CCD
            for c in i_f['ccd'].unique():
                b_ind = i_bi.loc[i_bi['ccd'] == c]
                f_ind = i_f.loc[(i_f['ccd'] == c) & (i_f['band'] == bx)]
                print(b_ind)
                print(f_ind)
                print(b_ind.index)
                print(f_ind.index)
                '''
                b_x1 = 
                b_x2 =
                f_x1 = 
                f_x2 = 
                pass
                # Data
                x1 = hdu[ext][a1:a2 + 1, a3:a4 + 1]
                x2 = hdu[ext][b1:b2 + 1, b3:b4 + 1]
                '''
        print('Exiting...')
    """

if __name__ == '__main__':
    # Directories where to look for biases and flats
    p_dir = '/archive_data/desarchive/'
    df_bias = pd.read_csv('set01_bias.csv')
    df_bias.columns = df_bias.columns.map(str.lower)
    df_bias['unitname'] = df_bias['unitname'].map(int)
    df_bias['path'] = map(lambda x: os.path.join(p_dir, x), df_bias['path'])
    #
    df_dflat = pd.read_csv('set01_dflat.csv')
    df_dflat.columns = df_dflat.columns.map(str.lower)
    df_dflat['unitname'] = df_dflat['unitname'].map(int)
    df_dflat['path'] = map(lambda x: os.path.join(p_dir, x), df_dflat['path'])

    # Call the methods
    t0 = time.time()
    method01(df_bias, df_dflat, fill_info=False)
    t1 = time.time()
    print('Elapsed time: {0:.2f}'.format((t1 - t0) / 60.))
