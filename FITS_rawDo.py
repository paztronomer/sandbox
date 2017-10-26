""" Script to perform statistics or other operations on raw exposures, using 
the DATASEC keyword (or others) to decide the sectioning from the whole focal
plane
"""

import socket
import argparse
import logging
import uuid
import copy
import numpy as np
import pandas as pd
import multiprocessing as mp
import fitsio

class FitsOpen():
    def __init__(self, fpath):
        f = fitsio.FITS(fpath)
        for x in f:
            print x
        self.ccd = np.copy(f.read())
        print self.ccd
        f.close()
    
    def get_ext():
        pass

class FitsImage():
    def __init__(self, fname=None, ftab=None):
        if not (fname is None):
            self.tab = np.array([fname])
        if not (ftab is None):
            self.tab = np.loadtxt(ftab, dtype=str, unpack=True)
        # In the future, change this dictionary for something more elegant
        # fpath, keys, outnm = kdic["fpath"], kdic["keys"], kdic["outnm"]
        # f = fitsio.FITS(fpath)
        # x_header = fitsio.read_header(fpath)
        # x_head = f[0].read_header()
        # x_hdu = f[0].read()
        # self.ccd = np.copy(x_hdu)
        # self.header = x_head
        # H = [x.read_header() for x in f]
        # f.close()
    
    def iter_f(self):
        for f in self.tab:
            FO = FitsOpen(f)
            FO.get_ext()

    def play_multiproc(self):
        proc_lst = []
        for i in range(10):
            proc = mp.Process(target=self.toy, args=(2,))
            proc_lst.append(proc)
        [x.start() for x in proc_lst]
        [x.join() for x in proc_lst]
    
    def toy(self, x):
        a = np.log(x*np.pi*1E45)
    
    
if __name__ == "__main__":
    # If file output, then use filename="some.txt" as argument
    logging.basicConfig(level=logging.DEBUG)
    logging.info(socket.gethostname())
    
    # Parse arguments
    t1 = "Perform some basic statistics on selected CCDs and save results"
    t1 += " on a csv file"
    par = argparse.ArgumentParser(description=t1)
    g1 = par.add_mutually_exclusive_group()
    h1 = "FITS file on which to operate"
    g1.add_argument("--fits", metavar="", help=h1)
    h2 = "List of paths to FITS files on which to operate"
    g1.add_argument("--path", metavar="", help=h2)
      
    nn = par.parse_args()
    d = dict()
    d["fname"] = nn.fits
    d["ftab"] = nn.path
    
    F = FitsImage(**d)
    F.iter_f()

    # Set multi-workers, process or what I need
    # d = mp.Manager().dict()
    # d["fpath"] = nn.fits
    # d["keys"] = nn.keyword
    # d["outnm"] = nn.o
    # proc = mp.Process(target=FitsImage, args=(d,))
    # proc.start()
    # proc.join()

