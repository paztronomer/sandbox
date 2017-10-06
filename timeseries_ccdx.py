"""Script to construct timeseries from a set of CCD images. It takes
different exposures, and given a CCD number, divide the area in rectangular
boxes for which statistics is calculated (under some normalization to allow
comparison among exposures).
Then, construct a 3D array (structured) and with this the time series is 
done.

Francisco Paz-Chinchon
DES/NCSA
"""

import os
import time
import socket
import argparse
import logging
import gc
import numpy as np
import pandas as pd
import scipy.signal as signal
import scipy.stats as stats
import fitsio
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
#import numpy.lib.recfunctions as recfunc


class FitsImage():
    def __init__(self,fpath):
        f = fitsio.FITS(fpath)
        #x_header = fitsio.read_header(fpath)
        x_head = f[0].read_header()
        x_hdu = f[0].read()
        self.ccd = np.copy(x_hdu)
        self.header = x_head
        f.close()
        #do I need to close the fits? M_hdu.close()


class Stack():
    def __init__(
            self, suffix=None, opt=None, table=None, loca=None,
            ccdnum=None, width_dim0=None, width_dim1=None, raw=None,
            region_x=None, region_y=None):
        if table is not None:
            self.df = pd.read_csv(table, sep=",")
        else:
            self.df = None
        if suffix is None:
            logging.info("\t-Suffix was not provided. Using \'d0xd1_PID\'")
            s = str(width_dim0) + "x" + str(width_dim1) + "_"
            s += str(os.getpid())
            self.suff = s
        else:
            self.suff = suffix 
        self.loca = loca
        self.opt = opt
        self.ccdnum = ccdnum
        self.w0 = width_dim0
        self.w1 = width_dim1
        self.raw = raw
        self.region_x = region_x
        self.region_y = region_y

    def norm_area(self, arr):
        """ Method to restrict (or not) the normalization value calculation 
        to a specific region on the CCD array. If using a raw CCD, must be
        cautious about it
        Inputs:
        - the whole CCD array
        Returns:
        - a subsection of the CCD array to be used for norm calculation
        """
        print arr.shape
        exit()
        if (self.raw):
            msg_w = "RAW image in use, with dimensions {0}".format(arr.shape)
            logging.warning(msg_w)
        if ( ~(self.region_x is None) and ~(self.region_y is None)):
            pass
        elif ( ~(self.region_x is None) and (self.region_y is None)):
            pass
        elif ( (self.region_x is None) and ~(self.region_y is None)):
            pass
        else:
            pass

    def one(self, root="/archive_data/desarchive"):
        """Runs the code over crosstalked CCDs, using paths to their 
        locations given in self.table
        """
        counter = 0
        for idx, row in self.df.iterrows():
            if np.equal(row["ccdnum"], self.ccdnum):
                gc.collect()
                aux = os.path.join(root, row["path"])
                aux = os.path.join(aux, row["filename"])
                fp = FitsImage(aux)
                obs = []
                obs += [int(fp.header["NITE"]), int(fp.header["EXPNUM"])]
                obs += [fp.header["BAND"].strip(), int(fp.header["CCDNUM"])]
                obs += [float(fp.header["EXPTIME"])]
                dt = np.dtype([("nite", "i4"), ("expnum", "i4"),
                               ("band", "|S10"), ("ccdnum", "i4"),
                               ("exptime", "f4")])
                obs = np.array([tuple(obs)], dtype=dt)
                if self.raw:
                    f1 = lambda x: int(x) - 1
                    dsec = fp.header["DATASEC"].strip().strip("[").strip("]")
                    dsec = map(f1, dsec.replace(":", ",").split(","))
                    # check if reads in reverse way
                    if (dsec[2] > dsec[3]) and (dsec[0] < dsec[1]):
                        fp_ccd_aux = fp.ccd[dsec[3] : dsec[2]+1, 
                                            dsec[0] : dsec[1]+1] 
                    elif (dsec[2] > dsec[3]) and (dsec[0] > dsec[1]):
                        fp_ccd_aux = fp.ccd[dsec[3] : dsec[2]+1, 
                                            dsec[1] : dsec[0]+1] 
                    elif (dsec[2] < dsec[3]) and (dsec[0] > dsec[1]):
                        fp_ccd_aux = fp.ccd[dsec[2] : dsec[3]+1, 
                                            dsec[1] : dsec[0]+1] 
                    elif (dsec[2] < dsec[3]) and (dsec[1] > dsec[0]):
                        fp_ccd_aux = fp.ccd[dsec[2] : dsec[3]+1, 
                                            dsec[0] : dsec[1]+1] 
                    else:
                        print "ERROR in DATA sectioning indices"
                        exit(1)
                    fp_ccd_aux += 1.
                else:
                    fp_ccd_aux = fp.ccd
                #create a tuple of arrays for the different sections,
                #iterating in dim0 and inside in dim1
                qx = Stats().quad(fp_ccd_aux, w0=self.w0, w1=self.w1)
                #here call the stats methods
                #Usual normalization is norm=np.median(fp.ccd))
                use4norm = self.norm_area(fp_ccd_aux)
                if self.opt == 1:
                    norm_x = np.median(use4norm)
                elif self.opt == 2:
                    norm_x = np.mean(use4norm)
                elif self.opt == 3:
                    norm_x = None
                else:
                    logging.error("Must select one of the available norm opt")
                    exit(1)
                qst = Stats().fill_it(qx, norm=norm_x)
                if counter == 0:
                    q3 = qst
                    obs3 = obs
                q3 = np.vstack((q3, qst))
                obs3 = np.vstack((obs3, obs))
                print counter
                counter += 1
        try:
            np.save("stat_{0:02}_{1}.npy".format(self.ccdnum, self.suff), q3)
            np.save("info_{0:02}_{1}.npy".format(self.ccdnum, self.suff), obs3)
            print "Done stats"
        except:
            logging.error("No computation was made. Exiting")
            exit(1)
        return True
    
    def two(self):
        """Run over images inside a folder
        """
        #walk in one depth level 
        DEPTH = 0
        c = 0
        for root, dirs, files in os.walk(self.loca):
            if root.count(os.sep) >= DEPTH:
                del dirs[:]
            for index, fits in enumerate(files):
                gc.collect()
                aux_f = os.path.join(self.loca, fits)
                f = FitsImage(aux_f) 
                if (self.ccdnum == f.header["CCDNUM"]):
                    #as nite is written by pipeline, must be constructed
                    nite = f.header["DATE-OBS"].strip().replace("-", "")[:8]
                    band = f.header["FILTER"].strip()[0]
                    h = [int(nite),f.header["EXPNUM"]]
                    h += [band,f.header["CCDNUM"]]
                    h += [f.header["EXPTIME"]]
                    dt = np.dtype([("nite", "i4"), ("expnum", "i4"),
                                   ("band", "|S10"), ("ccdnum", "i4"),
                                   ("exptime", "f4")])
                    hdr = np.array([tuple(h)], dtype=dt)
                    if self.raw:
                        f1 = lambda x: int(x) - 1
                        dsec = fp.header["DATASEC"].strip().strip("[").strip("]")
                        dsec = map(f1, dsec.replace(":", ",").split(","))
                        # check if reads in reverse way
                        if (dsec[2] > dsec[3]) and (dsec[0] < dsec[1]):
                            fp_ccd_aux = fp.ccd[dsec[3] : dsec[2]+1, 
                                                dsec[0] : dsec[1]+1] 
                        elif (dsec[2] > dsec[3]) and (dsec[0] > dsec[1]):
                            fp_ccd_aux = fp.ccd[dsec[3] : dsec[2]+1, 
                                                dsec[1] : dsec[0]+1] 
                        elif (dsec[2] < dsec[3]) and (dsec[0] > dsec[1]):
                            fp_ccd_aux = fp.ccd[dsec[2] : dsec[3]+1, 
                                                dsec[1] : dsec[0]+1] 
                        elif (dsec[2] < dsec[3]) and (dsec[1] > dsec[0]):
                            fp_ccd_aux = fp.ccd[dsec[2] : dsec[3]+1, 
                                                dsec[0] : dsec[1]+1] 
                        else:
                            print "ERROR in DATA sectioning indices"
                            exit(1)
                        fp_ccd_aux += 1.
                    else:
                        fp_ccd_aux = fp.ccd
                    qx = Stats().quad(f_ccd_aux, w0=self.w0, w1=self.w1)
                    use4norm = self.norm_area(fp_ccd_aux)
                    if self.opt == 1:
                        norm_x = np.median(use4norm)
                    elif self.opt == 2:
                        norm_x = np.mean(use4norm)
                    elif self.opt == 3:
                        norm_x = None
                    else:
                        logging.error("Must select a valid normalization")
                        exit(1)
                    qst = Stats().fill_it(qx, norm=norm_x)
                    if c == 0:
                        q3 = qst
                        hdr3 = hdr
                    q3 = np.vstack((q3, qst))
                    hdr3 = np.vstack((hdr3, hdr))
                    c += 1
        try:
            np.save("stat_{0:02}_{1}.npy".format(self.ccdnum, self.suff), q3)
            np.save("info_{0:02}_{1}.npy".format(self.ccdnum, self.suff), hdr3)
            print "Stats were successfully performed"
        except:
            logging.error("No computation was made. Exiting")
            exit(1)
        return True
                

class Stats():
    def quad(self, arr, w0=512, w1=1024):
        """Subdivides the image in smaller regions and returns a tuple
        of 2D arrays
        w0: width of dimension 0, for the subregion
        w1: width of dimension 1, for the subregion
        Also creates an ascii of the coordinates as a way to easily check
        positions
        """
        #limits for the subsections
        #the ending point is not in the array
        lim0 = np.arange(0, arr.shape[0], w0)
        lim1 = np.arange(0, arr.shape[1], w1)
        ss = [] #list of arrays
        cnt = 0
        with open("coord_{0}x{1}.csv".format(w0, w1), "w+") as f:
            m = "{0:<8}{1:<8}{2:<8}".format("index", "y_ini", "y_end")
            m += "{0:<8}{1:<8}\n".format("x_ini", "x_end")
            f.write(m)
            for j in lim0:
                for k in lim1:
                    ss.append(np.copy(arr[j : j+w0, k : k+w1]))
                    l = "{0:<8}{1:<8}{2:<8}".format(cnt, j, j + w0)
                    l += "{0:<8}{1:<8}\n".format(k, k + w1)
                    f.write(l)
                    cnt += 1
        return tuple(ss)

    def rms(self, arr):
        """returns RMS for ndarray
        """
        outrms = np.sqrt(np.mean(np.square(arr.ravel())))
        return outrms

    def uncert(self, arr):
        """calculates the uncertain in a parameter, as usually used in
        physics
        """
        ux = np.sqrt(np.mean(np.square(arr.ravel())) +
                     np.square(np.mean(arr.ravel())))
        return ux

    def mad(self, arr):
        return np.median(np.abs(arr - np.median(arr)))

    def entropy(self, arr):
        return stats.entropy(arr.ravel())

    def corr_random(self, arr):
        """correlate data with random 2D array
        """
        auxrdm = np.random.rand(arr.shape[0], arr.shape[1])
        auxrdm = auxrdm / np.mean(auxrdm)
        corr2d = signal.correlate2d(data, auxrdm, mode="same", boundary="symm")
        return corr2d

    def gaussian_filter(self, arr, sigma=1.):
        '''performs Gaussian kernel on image
        '''
        return scipy.ndimage.gaussian_filter(arr, sigma=sigma)

    def fill_it(self, tpl, norm=None):
        """For each one of the arrays in the tuple, perform statistics
        If no normalization value is given then result is not normalized
        Here the normalization is assumed as the division by a value
        """
        res = []
        dt = np.dtype([("norm", "f4"), ("med", "f4"),
                       ("avg", "f4"), ("med_n", "f4"),
                       ("avg_n", "f4"), ("rms_n", "f4"),
                       ("unc_n", "f4"), ("mad_n", "f4")])
        for x in tpl:
            if not isinstance(x, np.ndarray):
                logging.error("Not an array")
                exit(1)
            tmp = []
            if norm is None:
                norm = 1.
            tmp += [norm]
            tmp += [np.median(x), np.mean(x)]
            tmp += [np.median(x / norm), np.mean(x / norm)]
            tmp += [Stats().rms(x / norm), Stats().uncert(x / norm)]
            tmp += [Stats().mad(x / norm)]
            res.append(tuple(tmp))
        out = np.array(res, dtype=dt)
        return out


if __name__=="__main__":
    print socket.gethostname()
    #only temporary
    # tmp = "/work/devel/fpazch/calib_space/xtalkNoOversc_specter_y4e1/"
    # tmp += "CCD_1-2-3-6_noOversc"
    # tmp = "/work/devel/fpazch/calib_space/xtalk_specter_y4e1/CCD_2_3"
    tmp = "/work/devel/fpazch/calib_space/xtalkNoOversc_specter_y4e1/"
    tmp += "CCD_2-3_xtalked_ccd02NoOsc"
    #parser
    ecl = argparse.ArgumentParser(description="Time Series constructor")
    ecl.add_argument("-ccd", help="CCD number on which operate", metavar="",
                     type=int, default=3)
    ecl.add_argument("-norm", help="Normalization (1:med,2:avg,3:none)",
                     choices=[1,2,3], type=int)
    ecl.add_argument("--raw", help="Use if raw image with overscan",
                     action="store_true")
    g = ecl.add_mutually_exclusive_group()
    g.add_argument("--csv", help="Table with DB info (if needed)", metavar="")
    g.add_argument("--loc", help="Path to the CCD fits (if needed)",
                   metavar="", nargs="?", default=tmp, type=str)
    ecl.add_argument("--d0", help="Width of sub-boxes for dim0 (longer axis)",
                     metavar="", default=16, type=int)
    ecl.add_argument("--d1", help="Width of sub-boxes for dim1 (shorter axis)",
                     metavar="", default=128, type=int)
    ecl.add_argument("--suffix", help="Suffix for the output filenames",
                     metavar="")
    h_subx = "Shorter axis region to be used for calculate normalization value" 
    ecl.add_argument("--subx", help=h_subx, metavar="", nargs=2, type=int)
    h_suby = "Larger axis region to be used for calculate normalization value" 
    ecl.add_argument("--suby", help=h_suby, metavar="", nargs=2, type=int)
    nmsp = ecl.parse_args()
    # 
    print "\nInputs to the script:"
    v = vars(nmsp)
    for kx, vx in v.iteritems():
        print "\t{0:8}: {1}".format(kx, vx)
    #For skytemplates y4e1 I used: --csv redpixcor.csv
    kwin = {"table": nmsp.csv, "suffix": nmsp.suffix, "opt": nmsp.norm}
    kwin.update({"loca": nmsp.loc, "ccdnum": nmsp.ccd})
    kwin.update({"width_dim0": nmsp.d0, "width_dim1": nmsp.d1})
    kwin.update({"raw": nmsp.raw})
    kwin.update({"region_x": nmsp.subx, "region_y": nmsp.suby})
    if nmsp.csv is not None:
        Stack(**kwin).one()
    elif nmsp.loc is not None:
        Stack(**kwin).two()
