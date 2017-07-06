"""Script to construct timeseries from a set of CCD images
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
    def __init__(self,suffix=None,opt=None,table=None,loca=None,ccdnum=None,
                width_dim0=None,width_dim1=None):
        if table is not None:
            self.df = pd.read_csv(table,sep=",")
        else:
            self.df = None
        if suffix is None:
            logging.info("\t-Suffix was not provided. Using \'d0xd1_PID\'")
            s = str(width_dim0)+"x"+str(width_dim1)+"_"+str(os.getpid())
            self.suff = s
        else:
            self.suff = suffix 
        self.loca = loca
        self.opt = opt
        self.ccdnum = ccdnum
        self.w0 = width_dim0
        self.w1 = width_dim1

    def dbinfo(self):
        """
        select expnum,t_eff,b_eff,skybrightness,fwhm_asec from firstcut_eval 
           where expnum=567855
        select expnum,skypc00,skypc01,skypc02,skypc03,skyrms,skyfrac 
           from skyfit_qa where expnum=567855
        """
        pass

    def one(self,root="/archive_data/desarchive"):
        """Runs the code over crosstalked CCDs, using paths to their 
        locations given in self.table
        """
        counter = 0
        for idx,row in self.df.iterrows():
            if np.equal(row["ccdnum"],self.ccdnum):
                gc.collect()
                aux = os.path.join(root,row["path"])
                aux = os.path.join(aux,row["filename"])
                fp = FitsImage(aux)
                obs = []
                obs += [int(fp.header["NITE"]),int(fp.header["EXPNUM"])]
                obs += [fp.header["BAND"].strip(),int(fp.header["CCDNUM"])]
                obs += [float(fp.header["EXPTIME"])]
                dt = np.dtype([("nite","i4"),("expnum","i4"),("band","|S10"),
                            ("ccdnum","i4"),("exptime","f4")])
                obs = np.array([tuple(obs)],dtype=dt)
                #create a tuple of arrays for the different sections,
                #iterating in dim0 and inside in dim1
                qx = Stats().quad(fp.ccd,w0=self.w0,w1=self.w1)
                #here call the stats methods
                #Usual normalization is norm=np.median(fp.ccd))
                if self.opt == 1:
                    norm_x = np.median(fp.ccd)
                elif self.opt == 2:
                    norm_x = np.mean(fp.ccd)
                elif self.opt == 3:
                    norm_x = None
                else:
                    logging.error("Must select one of the available norm opt")
                    exit(1)
                qst = Stats().fill_it(qx,norm=norm_x)
                if counter == 0:
                    q3 = qst
                    obs3 = obs
                q3 = np.vstack((q3,qst))
                obs3 = np.vstack((obs3,obs))
                print counter
                counter += 1
        np.save("stat_{0}.npy".format(self.suff),q3)
        np.save("info_{0}.npy".format(self.suff),obs3)
        print "Done stats"
        return True
    
    def two(self):
        """Run over images inside a folder
        """
        #walk in one depth level 
        DEPTH = 0
        c = 0
        for root,dirs,files in os.walk(self.loca):
            if root.count(os.sep) >= DEPTH:
                del dirs[:]
            for index,fits in enumerate(files):
                gc.collect()
                aux_f = os.path.join(self.loca,fits)
                f = FitsImage(aux_f) 
                if (self.ccdnum == f.header["CCDNUM"]):
                    #as nite is written by pipeline, must be constructed
                    nite = f.header["DATE-OBS"].strip().replace("-","")[:8]
                    band = f.header["FILTER"].strip()[0]
                    h = [int(nite),f.header["EXPNUM"]]
                    h += [band,f.header["CCDNUM"]]
                    h += [f.header["EXPTIME"]]
                    dt = np.dtype([("nite","i4"),("expnum","i4"),
                                ("band","|S10"),("ccdnum","i4"),
                                ("exptime","f4")])
                    hdr = np.array([tuple(h)],dtype=dt)
                    qx = Stats().quad(f.ccd,w0=self.w0,w1=self.w1)
                    if self.opt == 1:
                        norm_x = np.median(f.ccd)
                    elif self.opt == 2:
                        norm_x = np.mean(f.ccd)
                    elif self.opt == 3:
                        norm_x = None
                    else:
                        logging.error("Must select a valid normalization")
                        exit(1)
                    qst = Stats().fill_it(qx,norm=norm_x)
                    if c == 0:
                        q3 = qst
                        hdr3 = hdr
                    q3 = np.vstack((q3,qst))
                    hdr3 = np.vstack((hdr3,hdr))
                    c += 1
        np.save("stat_{0}.npy".format(self.suff),q3)
        np.save("info_{0}.npy".format(self.suff),hdr3)
        print "Stats were successfully performed"
        return True
                

class Stats():
    def quad(self,arr,w0=512,w1=1024):
        """Subdivides the image in smaller regions and returns a tuple
        of 2D arrays
        w0: width of dimension 0, for the subregion
        w1: width of dimension 1, for the subregion
        Also creates an ascii of the coordinates as a way to easuly check
        positions
        """
        #limits for the subsections
        #the ending point is not in the array
        lim0 = np.arange(0,arr.shape[0],w0)
        lim1 = np.arange(0,arr.shape[1],w1)
        ss = [] #list of arrays
        cnt = 0
        with open("coord_{0}x{1}.csv".format(w0,w1),"w+") as f:
            m = "{0:<8}{1:<8}{2:<8}".format("index","y_ini","y_end")
            m += "{0:<8}{1:<8}\n".format("x_ini","x_end")
            f.write(m)
            for j in lim0:
                for k in lim1:
                    ss.append(np.copy(arr[j:j+w0,k:k+w1]))
                    l = "{0:<8}{1:<8}{2:<8}".format(cnt,j,j+w0)
                    l += "{0:<8}{1:<8}\n".format(k,k+w1)
                    f.write(l)
                    cnt += 1
        return tuple(ss)

    def rms(self,arr):
        """returns RMS for ndarray
        """
        outrms = np.sqrt(np.mean(np.square(arr.ravel())))
        return outrms

    def uncert(self,arr):
        """calculates the uncertain in a parameter, as usually used in
        physics
        """
        ux = np.sqrt(np.mean(np.square(arr.ravel())) +
                    np.square(np.mean(arr.ravel())))
        return ux

    def mad(self,arr):
        return np.median(np.abs(arr - np.median(arr)))

    def entropy(self,arr):
        return stats.entropy(arr.ravel())

    def corr_random(self,arr):
        """correlate data with random 2D array
        """
        auxrdm = np.random.rand(arr.shape[0],arr.shape[1])
        auxrdm = auxrdm/np.mean(auxrdm)
        corr2d = signal.correlate2d(data,auxrdm,mode="same",boundary="symm")
        return corr2d

    def gaussian_filter(self,arr,sigma=1.):
        '''performs Gaussian kernel on image
        '''
        return scipy.ndimage.gaussian_filter(arr,sigma=sigma)

    def fill_it(self,tpl,norm=None):
        """For each one of the arrays in the tuple, perform statistics
        If no normalization value is given then result is not normalized
        Here the normalization is assumed as the division by a value
        """
        res = []
        dt = np.dtype([("med","f4"),("avg","f4"),
                    ("med_n","f4"),("avg_n","f4"),
                    ("rms_n","f4"),("unc_n","f4"),
                    ("mad_n","f4")])
        for x in tpl:
            if not isinstance(x,np.ndarray):
                logging.error("Not an array")
                exit(1)
            tmp = []
            if norm is None:
                norm = 1.
            tmp += [np.median(x),np.mean(x)]
            tmp += [np.median(x/norm),np.mean(x/norm)]
            tmp += [Stats().rms(x/norm),Stats().uncert(x/norm)]
            tmp += [Stats().mad(x/norm)]
            res.append(tuple(tmp))
        out = np.array(res,dtype=dt)
        return out


if __name__=="__main__":
    print socket.gethostname()
    #only temporary
    tmp = "/work/devel/fpazch/calib_space/xtalkNoOversc_specter_y4e1/"
    tmp += "CCD_1-2-3-6_noOversc"
    #parser
    ecl = argparse.ArgumentParser(description="Time Series constructor")
    ecl.add_argument("-norm",help="Normalization (1:med,2:avg,3:none)",
                    choices=[1,2,3],type=int)
    g = ecl.add_mutually_exclusive_group()
    g.add_argument("--csv",help="Table with DB info (if needed)",metavar="")
    g.add_argument("--loc",help="Path to the CCD fits (if needed)",metavar="",
                default=tmp)
    ecl.add_argument("--ccdnum",help="CCD number on which operate",metavar="",
                    type=int,default=3)
    ecl.add_argument("--d0",help="Width of sub-boxes for dim0 (longer axis)",
                    metavar="",default=16,type=int)
    ecl.add_argument("--d1",help="Width of sub-boxes for dim1 (shorter axis)",
                    metavar="",default=128,type=int)
    ecl.add_argument("--suffix",help="Suffix for the output filenames",
                    metavar="")
    nmsp = ecl.parse_args()
    #For skytemplates y4e1 I used: --csv redpixcor.csv
    kwin = {"table":nmsp.csv,"suffix":nmsp.suffix,"opt":nmsp.norm}
    kwin.update({"loca":nmsp.loc,"ccdnum":nmsp.ccdnum})
    kwin.update({"width_dim0":nmsp.d0,"width_dim1":nmsp.d1})
    if nmsp.csv is not None:
        Stack(**kwin).one()
    elif nmsp.loc is not None:
        Stack(**kwin).two()
