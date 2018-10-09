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
import argparse
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
    fnm_f1, fnm_f2, fnm_b1, fnm_b2, band, ccd, nite, box_side, ext, w_amp, h_amp = inlist 
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
             outname=None,
             gainout=None,
             box_side=256,
             ext=0,
             w_amp=1024,
             h_amp=2048,):
    """ Calculates gain changes  
    Inputs
    - tab: dataframe, must contain columns path, reqnum, unitname (=nite)
    """
    if (outname is not None):
        fill_info = True
    # Define pairs of bias, domeflats. Use bias simply by exposure number 
    # order. 
    # Do calculation by band, by CCD, by amplifier
    # Params to read FITS section
    a1, a2 = (h_amp - 1) - int(box_side / 2), (h_amp - 1) + int(box_side / 2)
    a3 = (int(w_amp / 2) - 1) - int(box_side / 2), 
    a4 = (int(w_amp / 2) - 1) + int(box_side / 2) 
    b1, b2 = (h_amp - 1) - int(box_side / 2), (h_amp - 1) + int(box_side / 2)
    aux_2amp_w = (w_amp * 2) - int(w_amp / 2) - 1
    b3 = aux_2amp_w - int(box_side / 2)
    b4 = aux_2amp_w + int(box_side / 2) 
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
                        box_side,
                        ext,
                        w_amp,
                        h_amp,
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
    print('Writen {0}'.format(gainout))
    return True     

if __name__ == '__main__':
    txt_desc = 'Gain calculation using the Janesick method. Note we suppose'
    txt_desc += ' each CCD is DES-like, two amplifiers side by side'
    txt_desc += ' attached by the longer axis.'
    txt_desc += ' For method see:'
    txt_desc += ' http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?findgain'
    arg = argparse.ArgumentParser(description=txt_desc)
    #
    h0 = 'Table containing the BIAS paths to the xtalked CCDs. Columns: PATH,'
    h0 += 'REQNUM,UNITNAME. Format: CSV'
    arg.add_argument('bias', help=h0)
    h1 = 'Table containing the DOME FLAT paths to the xtalked CCDs. Columns: '
    h1 += 'PATH,REQNUM,UNITNAME. Format: CSV'
    arg.add_argument('dflat', help=h1)
    default_prefix = '/archive_data/desarchive'
    h2 = 'Prefix to be attached before every BIAS/DOME FLAT path. For no'
    h2 += ' prefix use the variable followed by blank space.'
    h2 += ' Default: {0}'.format(default_prefix)
    arg.add_argument('--pre', help=h2, default=default_prefix)
    h3 = 'Output name for the file containing basic info and gain from each'
    h3 += ' image header. If no name is given, no file will be generated'
    arg.add_argument('--head', help=h3)
    default_gain_fnm = 'gain_PID{0}.csv'.format(os.getpid())
    h4 = 'Output name for the gain calculation CSV table. Default:'
    h4 += ' {0}'.format(default_gain_fnm)
    arg.add_argument('--res', help=h4, default=default_gain_fnm)
    # 
    gain_box = 256 
    ext = 0
    width_amp = 1024
    height_amp = 2048
    h5 = 'Box side in pix for GAIN measurement. Default: {0}'.format(gain_box)
    arg.add_argument('--box', help=h5, default=gain_box)
    h6 = 'Extension from where to read the image. Default: {0}'.format(ext)
    arg.add_argument('--ext', help=h6, default=ext)
    h7 = 'Dimensions of each amplifier. Format: long-axis short-axis'
    h7 += ' Default: {0} {1}'.format(height_amp, width_amp)
    arg.add_argument('--dim', help=h7, default=[height_amp, width_amp], 
                     nargs=2)
    #
    arg = arg.parse_args()
    #
    # Before to run, prepare the dataframes
    #
    """
    outname='header_info_tab.csv',
    gainout='calc_gain_perAmp.csv'
    box_side=256 
    ext=0
    w_amp=1024
    h_amp=2048
    """
    # Directories where to look for biases and flats
    df_bias = pd.read_csv(arg.bias)
    df_bias.columns = df_bias.columns.map(str.lower)
    df_bias['unitname'] = df_bias['unitname'].map(int)
    df_bias['path'] = map(lambda x: os.path.join(arg.pre, x), 
                          df_bias['path'])
    #
    df_dflat = pd.read_csv(arg.dflat)
    df_dflat.columns = df_dflat.columns.map(str.lower)
    df_dflat['unitname'] = df_dflat['unitname'].map(int)
    df_dflat['path'] = map(lambda x: os.path.join(arg.pre, x), 
                           df_dflat['path'])
    # Call the methods
    t0 = time.time()
    method01(df_bias, df_dflat, outname=arg.head, gainout=arg.res,
             box_side=arg.box, ext=arg.ext, 
             h_amp=arg.dim[0], w_amp=arg.dim[1],)
    t1 = time.time()
    print('Elapsed time: {0:.2f}'.format((t1 - t0) / 60.))
