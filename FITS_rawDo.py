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
        """ Open the FITS and iterate over within extensions. Note once close()
        is called, data cannot be accessed 
        """
        self.ccd = fitsio.FITS(fpath)
        logging.info("Working on FITS file: {0}".format(fpath))
    
    def get_ext(self, npix=120):
        """ Go through the different extensions, select those not being FOCUS
        CCDs and, after performing masking around left and right edges (tape
        bumps), get some basic statistics
        """
        expnum = self.ccd[0].read_header()["EXPNUM"]
        band = self.ccd[0].read_header()["FILTER"].split()[0]
        # Get which header extensions are science. Remember first extension
        # has not EXTNAME keyword
        extnm, aux_idx = [], []
        for idx, x in enumerate(self.ccd):
            try:
                extnm.append(x.read_header()["EXTNAME"].strip())
                aux_idx.append(idx)
            except:
                pass
        extnm, aux_idx = np.array(extnm), np.array(aux_idx)
        condition = np.flatnonzero(np.core.defchararray.find(extnm, "F"))
        extnm = extnm[condition]
        aux_idx = aux_idx[condition]
        # Using the indices of the extensions not being FOCUS CCDs, perform the
        # simple statistics
        res = []
        for xtens in aux_idx:
            tmp_h = self.ccd[xtens].read_header()
            tmp_data = self.ccd[xtens].read()
            # Double checking
            if not (tmp_h["EXTNAME"].strip() in extnm):
                logging.error("EXTNAME not in list")
                exit(1)
            # Use DATASEC{A,B} to get order and to trim. Save order
            dsA = tmp_h["DATASECA"].strip("[").strip("]").replace(":", ",")
            dsB = tmp_h["DATASECB"].strip("[").strip("]").replace(":", ",")
            dsA = map(int, dsA.split(","))
            dsB = map(int, dsB.split(","))
            # Order AB = 1, BA = -1
            if (dsA[0] < dsB[0]):
                order = 1
                dsA[0] += npix
                dsB[1] -= npix
            else:
                order = -1
                dsB[0] += npix
                dsA[1] -= npix
            ampA = tmp_data[dsA[2] : dsA[3], dsA[0] : dsA[1]]
            ampB = tmp_data[dsB[2] : dsB[3], dsB[0] : dsB[1]]
            # Basic statistics
            stA = [expnum, band, tmp_h["CCDNUM"], "A", order] 
            stA += self.basic_stat(ampA) 
            stB = [expnum, band, tmp_h["CCDNUM"], "B", order] 
            stB += self.basic_stat(ampB)
            # Save results
            res.append(stA)
            res.append(stB)
            # Order or reading
            # NOTE: if need to read in reverse order, better than [::-1] is to
            #       use np.flip(arr, axis)
        # Call the FITS instace closing
        self.close_fits()
        logging.info("Done with calculations onf expnum: {0}".format(expnum))
        return res

    def basic_stat(self, x):
        """ Method to perform basic statistics on the array
        """
        s1 = np.median(x.ravel())
        s2 = np.mean(x.ravel())
        s3 = np.sqrt(np.mean(np.square(x.ravel()))) 
        s4 = np.sqrt(np.mean( np.square(x.ravel()) ) + 
                     np.square( np.mean(x.ravel()) )) 
        s5 = np.var(x)
        s6 = np.min(x.ravel())
        s7 = np.max(x.ravel())
        s8 = np.percentile(x.ravel(), 5)
        s9 = np.percentile(x.ravel(), 25)
        s10 = np.percentile(x.ravel(), 75)
        s11 = np.percentile(x.ravel(), 95)
        return [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11]

    def close_fits(self):
        """  Closes the open fitsio.FITS instance
        """
        self.ccd.close()
        return True


class Loader():
    def __init__(self, fname=None, ftab=None, suffix=None):
        """ Loader of the files or tables, to be passed to the class who 
        actually do the job
        """
        if not (fname is None):
            self.tab = np.array([fname])
        if not (ftab is None):
            self.tab = np.loadtxt(ftab, dtype=str, unpack=True)
        if (suffix is None):
            self.suffix = str(uuid.uuid4())
        else:
            self.suffix=suffix

    def iter_f(self):
        """ Simple iterator to call the FITS loader and save the results  
        """
        tmp = []
        for f in self.tab:
            FO = FitsOpen(f)
            tmp += FO.get_ext()
        # The sum of allresults must be saved in CSV
        logging.info("Checkpoint")
        np.save("ccdstatistics_checkpoint.npy", np.array(tmp))
        aux_tmp = zip(*tmp)
        cols = ["expnum", "band", "ccdnum", "amp", "amp_order", 
                "med", "avg", "rms", "uncert", "var", "min", "max", 
                "p5", "p25", "p75", "p95"]
        d = dict()
        for idx, item in enumerate(cols):
            d[item] = aux_tmp[idx]
        outnm = "ccdStat_{0}.csv".format(self.suffix)
        pd.DataFrame(data=d).to_csv(outnm, index=False, header=True)
        logging.info("Results were saved: {0}".format(outnm))
        return True

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
    
    L = Loader(**d)
    L.iter_f()

    # Set multi-workers, process or what I need
    # d = mp.Manager().dict()
    # d["fpath"] = nn.fits
    # d["keys"] = nn.keyword
    # d["outnm"] = nn.o
    # proc = mp.Process(target=FitsImage, args=(d,))
    # proc.start()
    # proc.join()

