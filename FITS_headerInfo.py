""" Code to get information from header, iterating over all FITS extensions
"""

import socket
import argparse
import logging
import uuid
import numpy as np
import multiprocessing as mp
import fitsio


class FitsImage():
    def __init__(self, fpath=None, keys=None, outnm=None):
        f = fitsio.FITS(fpath)
        # x_header = fitsio.read_header(fpath)
        x_head = f[0].read_header()
        x_hdu = f[0].read()
        self.ccd = np.copy(x_hdu)
        self.header = x_head
        if (outnm is None):
            self.out = str(uuid.uuid4()) + ".csv"
        else:
            if (outnm[-4:] == ".csv"):
                self.out = outnm
            else:
                self.out = outnm + ".csv"
        f.close()

    def get_kw(self):
        """ To get the keywords and contruct the output file
        """
        for i in range(100):
            print self.header

if __name__ == "__main__":
    # If file output, then use filename="some.txt" as argument
    logging.basicConfig(level=logging.DEBUG)
    logging.info(socket.gethostname())

    # Parse arguments
    t1 = "Code to get header keyword info from all FITS extensions and save"
    t1 += " into a file."
    eel = argparse.ArgumentParser(description=t1)
    t2 = "Filename of the FITS file on which to run"
    eel.add_argument("fits", help=t2)
    t3 = "List of space-separated values to be readed from Header of different"
    t3 += "FITS extensions"
    eel.add_argument("keyword", help=t3, nargs="?")
    t4 = "Output name for the generated CSV table of keywords"
    eel.add_argument("-o", help=t4, metavar="")
    # t4 = "Number of processes to be started (multiprocessing). Default: 4"
    # eel.add_argument("-p", help=t4, metavar="", default=4, type=int)
    nn = eel.parse_args
    d = dict()
    d["fpath"] = nn.fits
    d["keys"] = nn.keyword

    # Initialize

    # Set multi-workers, process or what I need
    proc = mp.Process(target=FitsImage, args=(kwargs=d))
    proc.start()
    proc.join()

