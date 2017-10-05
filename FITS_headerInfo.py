""" Code to get information from header, iterating over all FITS extensions
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


class FitsImage():
    def __init__(self, fpath=None, keys=None, outnm=None, science=None,
                 printX=None):
        # In the future, change this dictionary for something more elegant
        # fpath, keys, outnm = kdic["fpath"], kdic["keys"], kdic["outnm"]
        f = fitsio.FITS(fpath)
        # x_header = fitsio.read_header(fpath)
        # x_head = f[0].read_header()
        # x_hdu = f[0].read()
        # self.ccd = np.copy(x_hdu)
        # self.header = x_head
        H = [x.read_header() for x in f]
        self.header = H
        self.kw = keys
        if (outnm is None):
            self.out = str(uuid.uuid4()) + ".csv"
        else:
            if (outnm[-4:] == ".csv"):
                self.out = outnm
            else:
                self.out = outnm + ".csv"
        self.ccdsci = science
        self.printX = printX
        f.close()
    
    def get_kw(self, initial_ext=1):
        """ To get the keywords and contruct the output file
        """
        # If printing was selected
        if (not (self.printX is None)):
            for hx in self.header[initial_ext:]:
                if (hx["EXTNAME"].strip() == self.printX.upper()):
                    print "{0}\n{1}\n{2}\n".format("="*80, hx, "="*80)
                    check_print = True
                    exit(0)
            # If the code didn't exit at this point, then the extension 
            # wasn't found
            m_print = "Extension name: {0} was not found".format(self.printX)
            logging.error(m_print)
            exit(1)
        # For each of the extensions, construct the table
        # List of keywords to be extracted, different for Focus and for Science
        if (len(self.kw) == 0):
            # If no one key was required, then all keywords will be 
            # returned, assuming all headers has same entries
            cont = True
            aux_i = initial_ext
            while cont:
                j = self.header[aux_i]
                if ((self.ccdsci) and (j["EXTNAME"][0].strip() != "F")):
                    self.kw = [x for x in j]
                    cont = False
                elif ((not self.ccdsci) and (j["EXTNAME"][0].strip() == "F")):
                    self.kw = [z for z in j]
                    cont = False
                aux_i += 1
        # Fill a list with the values of the header, making the distinction if
        # Science or Focus was needed
        all_ok = True
        tmp_x = []
        r_space = lambda x: x.strip() if isinstance(x, str) else x
        for idx, ext in enumerate(self.header[initial_ext:]):
            if ((self.ccdsci) and (ext["EXTNAME"][0].strip() != "F")):
                try:
                    tab_x = [ext[key_i] for key_i in self.kw]
                    tab_x = map(r_space, tab_x)
                    tmp_x.append(tab_x)
                except Exception as e:
                    msg = "{0} / EXTENSION: {1}".format(e, ext["EXTNAME"])
                    logging.error(msg)
                    all_ok = False
            elif ((not self.ccdsci) and (ext["EXTNAME"][0].strip() == "F")):
                tab_x = [ext[key_i] for key_i in self.kw]
                tab_x = map(r_space, tab_x)
                tmp_x.append(tab_x)
        # Check if all extensions were successfully queried
        if (not all_ok):
            m_out = "="*12
            m_out += "Exiting. Not all extensions were succesfully read."
            m_out += "="*12
            logging.error(m_out)
            exit(1)
        # Construct the dataframe and save to file
        tmp_y = zip(*tmp_x)
        df = pd.DataFrame(data=dict(zip(self.kw, tmp_y)))
        logging.info("Saving CSV: {0}".format(self.out))
        df.to_csv(self.out, header=True, index=False)
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
    t1 = "Code to get header keyword info from all FITS extensions and save"
    t1 += " into a CSV file. Do not consider the first extension. Using DECam."
    eel = argparse.ArgumentParser(description=t1)
    t2 = "Filename of the FITS file on which to run"
    eel.add_argument("fits", help=t2)
    t3 = "List of space-separated values to be readed from header of different"
    t3 += " FITS extensions. Default: all the keywords"
    eel.add_argument("keyword", help=t3, nargs="*")
    t4 = "Output name for the generated CSV table of keywords"
    eel.add_argument("-o", help=t4, metavar="")
    t5 = "Flag to use Focus CCDs instead of Science."
    eel.add_argument("-f", help=t5, action="store_false")
    t6 = "Want to print an extension header? Input its \'EXTNAME\'. Script"
    t6 += " will exit after printing."
    eel.add_argument("-p", help=t6, metavar="")
    # t4 = "Number of processes to be started (multiprocessing). Default: 4"
    # eel.add_argument("-p", help=t4, metavar="", default=4, type=int)
    nn = eel.parse_args()
    # Initialize
    d = dict()
    d["fpath"] = nn.fits
    d["keys"] = nn.keyword
    d["outnm"] = nn.o
    d["printX"] = nn.p
    d["science"] = nn.f

    F = FitsImage(**d)
    F.get_kw()
    
    # Set multi-workers, process or what I need
    # d = mp.Manager().dict()
    # d["fpath"] = nn.fits
    # d["keys"] = nn.keyword
    # d["outnm"] = nn.o
    # proc = mp.Process(target=FitsImage, args=(d,))
    # proc.start()
    # proc.join()

