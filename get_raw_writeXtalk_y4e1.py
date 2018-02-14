"""To get paths for Y4E1 skytemplates source raw exposures and write 
xtalk decam script
"""

import os
import sys
import logging
import time 
import numpy as np
import pandas as pd
import argparse
import tempfile
import fitsio
import despydb.desdbi as desdbi

class Toolbox(object):    
    @classmethod
    def dbquery(cls, toquery, outdtype, dbsection='db-desoper', help_txt=False):
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
        outtab = np.rec.array(rows, dtype=zip(cols, outdtype))
        return outtab

    @classmethod
    def rawexp(cls, expnum, root="/archive_data/desarchive"):
        """With the expnum, returns the full path
        """
        q = "select fai.path,fai.filename,fai.compression"
        q += " from file_archive_info fai, exposure e"
        q += " where e.expnum={0}".format(expnum)
        q += " and e.filename=fai.filename"
        dt = ["a100","a50","a10"]
        res = Toolbox.dbquery(q, dt)
        if res["path"].shape[0] > 1:
            logging.error("More than one occurrence for {0}".format(expnum))
            exit(1)
        p = os.path.join(root,res["path"][0])
        p = os.path.join(p,res["filename"][0] + res["compression"][0])
        aux_fnm = res["filename"][0]
        aux_fnm = aux_fnm[:aux_fnm.find(".fits")]
        return p, aux_fnm

if __name__=="__main__":
    gral = "Script to use a set of exposure numbers (EXPNUM) to query the DB,"
    gral += " and with the available images construct a bash script to run"
    gral += " crosstalk. If the \'modified\' crosstalk version is employed"
    gral += " (RGruendl version), then CCD02 overscan can be turned off."
    gral += " \nIMPORTANT: header updater is NOT the last version."
    arg = argparse.ArgumentParser(description=gral) 
    h1 = "Exposure list, having \'EXPNUM\' as column name. File can contain"
    h1 += " mutiple columns"
    arg.add_argument("exposures", help=h1)
    h2 = "Flag to turn-off the overscan subtraction on CCD02 only"
    arg.add_argument("-n", help=h2, action="store_true")
    h3 = "Suffix to add at the end of the output filenames"
    arg.add_argument("-s", help=h3, metavar="")
    h4 = "Filename of the crosstalk coefficients matrix."
    h4 += " Default: DECam_20130606.xtalk"
    arg.add_argument("-c", help=h4, metavar="")
    h5 = "List of space-separated CCD numbers to be xtalked."
    h5 += " Default is to use all"
    arg.add_argument("-l", help=h5, metavar="", nargs="*", type=int)
    h6 = "Print a sample command line, with the selected parameters"
    arg.add_argument("-p", help=h6, action="store_true")
    h7 = "Suffix to add to the bash/csh script filename"
    arg.add_argument("-o", help=h7, metavar="")
    h8 = "Method to evaluate overscan, for M=0 the whole row is used,"
    h8 += " otherwise M represents the number of columns to be used."
    h8 += " M=[-50, -1] cubic spline; M=0 line by line; M=[1, 50] Legendre."
    h8 += " Default: 0. This option is equivalent to \'overscanfunction\'"
    arg.add_argument("-of", help=h8, metavar="", default=0, type=int)
    h9 = "Index to be used as the order for the Legendre polynomial, when"
    h9 += " selected on the \'overscanfunction\'. Order=[1, 6]. Default: 1."
    h9 += " This option is equivalent to \'overscanorder\'"
    arg.add_argument("-oo", help=h9, metavar="", default=1, type=int)
    val = arg.parse_args()
    # Arguments
    aux_null = val.n
    if isinstance(val.l, list):
        ccds = ",".join(map(str, val.l))
    else:
        ccds = val.l
    # Reading table of expnums 
    df = pd.read_csv(val.exposures)
    wr = []
    for idx,row in df.iterrows():
        route,aux = Toolbox.rawexp(row["EXPNUM"])
        one = "DECam_crosstalk"
        if (val.s is None):
            two = " " + aux + "_%02d.fits"
        else:
            two = " " + aux + "_%02d_{0}.fits".format(val.s)
        two += " -overscansample 1 -overscantrim 5"
        # adding overscanfunction
        two += " -overscanfunction {0}".format(val.of)
        # adding overscanorder
        two += " -overscanorder {0}".format(val.oo)
        if aux_null:
            # For no overscan CCD2 only
            two += " -null_overscan_ccd2"
        two += " -replace /archive_data/desarchive/OPS/config/20170531/"
        two += "finalcut/20170531_DES_header_update.20140303"
        two += " -verbose 3  -photflag 1"
        if (val.c is None):
            two += " -crosstalk /archive_data/desarchive/OPS/cal/xtalk/"
            two += "20130606/DECam_20130606.xtalk"
        else:
            two += " -crosstalk {0}".format(val.c)
        if (ccds is not None):
            two += " -ccdlist {0}".format(ccds)
        two += " > log/log_{0:08}_pid{1}.txt".format(row["EXPNUM"], os.getpid())
        out = one + " " + route + two
        if (idx == 0) and (val.p):
            print "\nSample command line:\n{0}\n".format(out)
        wr.append(out)
    #
    # Temporary file (and filename) generator. Alternative to os.getpid()
    aux_name = tempfile.NamedTemporaryFile(
        mode="w+b", 
        prefix="xtalk_", 
        dir=os.getcwd(),
        delete=True)
    # As we only need the name and not the file, set delete to True
    aux_name = aux_name.name
    if (val.o is None):
        aux_name += ".csh"
    else:
        aux_name += "_{0}.csh".format(val.o)
    #
    with open(aux_name, "wr") as w:
        w.write("#!/bin/csh \n")
        for x in wr:
            w.write(x)
            w.write("\n")
    print "Done!\n{0}".format(time.ctime())
