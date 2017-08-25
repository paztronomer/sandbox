"""To get paths for Y4E1 skytemplates source raw exposures and write 
xtalk decam script
"""

import os
import sys
import logging
import time 
import numpy as np
import pandas as pd
import fitsio
import despydb.desdbi as desdbi

class Toolbox(object):    
    @classmethod
    def dbquery(cls,toquery,outdtype,dbsection='db-desoper',help_txt=False):
        '''the personal setup file .desservices.ini must be pointed by desfile
        DB section by default will be desoper
        '''
        desfile = os.path.join(os.getenv('HOME'),'.desservices.ini')
        section = dbsection
        dbi = desdbi.DesDbi(desfile,section)
        if help_txt: help(dbi)
        cursor = dbi.cursor()
        cursor.execute(toquery)
        cols = [line[0].lower() for line in cursor.description]
        rows = cursor.fetchall()
        outtab = np.rec.array(rows,dtype=zip(cols,outdtype))
        return outtab

    @classmethod
    def rawexp(cls,expnum,root="/archive_data/desarchive"):
        """With the expnum, returns the full path
        """
        q = "select fai.path,fai.filename,fai.compression"
        q += " from file_archive_info fai, exposure e"
        q += " where e.expnum={0}".format(expnum)
        q += " and e.filename=fai.filename"
        dt = ["a100","a50","a10"]
        res = Toolbox.dbquery(q,dt)
        if res["path"].shape[0] > 1:
            logging.error("More than one occurrence for {0}".format(expnum))
            exit(1)
        p = os.path.join(root,res["path"][0])
        p = os.path.join(p,res["filename"][0]+res["compression"][0])
        aux_fnm = res["filename"][0]
        aux_fnm = aux_fnm[:aux_fnm.find(".fits")]
        return p,aux_fnm

if __name__=="__main__":
    fname = sys.argv[1]
    print "Input table of expnums: {0}".format(fname)
    df = pd.read_csv(fname)
    wr = []
    for idx,row in df.iterrows():
        route,aux = Toolbox.rawexp(row["EXPNUM"])
        one = "DECam_crosstalk"
        two = " " + aux + "%02d_xtalked.fits" 
        # Two options: no overscan and the traditional
        # two += + "_noOver"
        two += " -overscansample 1  -overscanfunction 0  -overscantrim 5"
        two += " -replace /archive_data/desarchive/OPS/config/20170531/"
        two += "finalcut/20170531_DES_header_update.20140303"
        two += " -verbose 3  -photflag 1"
        two += " -crosstalk /archive_data/desarchive/OPS/cal/xtalk/20130606/"
        two += "DECam_20130606.xtalk"
        two += " -ccdlist 2,3"
        two += " > log/log_{0:08}.txt".format(row["EXPNUM"])
        out = one + " " + route + two
        wr.append(out)
    with open("{0}.csh".format(os.getpid()),"wr") as w:
        w.write("#!/bin/csh \n")
        for x in wr:
            w.write(x)
            w.write("\n")
    print "Done!"
