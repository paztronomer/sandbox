""" Simple code to change the value of a keyword
"""


import os
import time
import logging
import argparse
import numpy as np
import multiprocessing as mp
import astropy.io.fits as apy_fits
import fitsio

logging.basicConfig(
                    #filename = 'out.log',
                    level = logging.DEBUG,
                    format = "%(asctime)s - %(levelname)s - %(message)s"
                   )
def astropy_mod(kw):
    """ Astropy function to update header. Fitsio does not seems to be so 
    direct
    """
    fnm, key_ls , val_ls, ext_ls = kw
    with apy_fits.open(fnm, 'update') as f:
        for idx, xt in enumerate(ext_ls):
            f[xt].header[key_ls[idx]] = val_ls[idx]
    return True

def mod_key(k_ls=None, str_ini=None, str_end=None, path_ls=None):
    """ Modifies the key
    Inputs
    - k_list: list of keywords to modify
    - str_ini: list of the initial keyword values
    - str_end: list of the values to be set in replacement
    - path_ls: list of parent paths were the fits files are located. Recursive
    search will be performed to reach all levels
    """
    # Iteratively search
    fnm_ls = []
    for p in path_ls:
        for root, dirlist, files in os.walk(p):
            for f in filter(lambda x: '.fits' in x, files):
                fnm_ls.append(os.path.join(root, f))
    # Save a list to then run in parallel the header modification
    aux_mod = []
    for f in fnm_ls:
        keyw_tmp, ext_tmp, val_tmp = [], [], []
        try:
            hdu = fitsio.FITS(f, 'rw')
            n_ext = len(hdu)
            for xt in range(n_ext):
                # Try reading the header
                try:
                    h = hdu[xt].read_header()
                    # Try get the desired  keywords
                    for idx, k in enumerate(k_ls):
                        try:
                            # If keyword is found, then change sub-string
                            aux = h[k]
                            aux = aux.replace(str_ini[idx], str_end[idx])
                            # Keywords, extensions
                            keyw_tmp.append(k)
                            ext_tmp.append(xt)
                            val_tmp.append(aux)
                        except:
                            logging.info('Key not found:{0} {1}'.format(k, f))
                except:
                    logging.error('Unable to read the header: {0}'.format(f))
            
            # Fill auxiliary list for header modify
            if ((len(keyw_tmp) > 0) and (len(val_tmp) > 0) 
                and (len(ext_tmp) > 0)):
                aux_mod.append([f, keyw_tmp, val_tmp, ext_tmp])
            # Closing HDU
            hdu.close()
        except:
            raise
    # Using astropy, change the files in parallel
    P1 = mp.Pool(processes=mp.cpu_count())
    res= P1.map(astropy_mod, aux_mod)
    P1.close()

    return True

if __name__ == '__main__':
    t_intro = 'Simple script to update keywords in a set of FITS files, in'
    t_intro += ' parallel'
    arg = argparse.ArgumentParser(description=t_intro)
    #
    t0 = 'Space separated list of keywords to be modified'
    arg.add_argument('-k', metavar='list', help=t0, nargs='+')
    t1 = 'Space separated list of strings to look for. One per input keyword.'
    t1 += ' If multiple words compose the string, use \' \' to enclose them.'
    t1 += ' Example: \'SDSS z filter c0003\'' 
    arg.add_argument('-s0', metavar='initial', help=t1, nargs='+')
    t2 = 'Space separated list of strings used as replacement.'
    t2 += ' One per input keyword'
    t2 += ' If multiple words compose the string, use \' \' to enclose them.'
    t2 += ' Example: \'WISE W1 filter\'' 
    arg.add_argument('-s1', metavar='final', help=t2, nargs='+')
    t3 = 'Space separated list of parent paths were to recursive loof for FITS'
    t3 += ' files. All keywords will be searched on every found file'
    arg.add_argument('-p', metavar='path', help=t3, nargs='+')
    #
    v = arg.parse_args()
    # Call the method
    mod_key(k_ls=v.k, str_ini=v.s0, str_end=v.s1, path_ls=v.p)


    #
    #k_ls = ['BAND', 'FILTER']
    #str_end = ['W', 'W']
    #str_ini = ['N964', 'N964']
    #path_ls = ['20170815t0222-r3572', 
    #           '20170815t0824-r3571', 
    #           '20170906t0222-r3573']
