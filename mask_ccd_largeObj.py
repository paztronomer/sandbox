""" Script to mask CCD by CCD, based on the outputs from SExtractor. The idea
behind can be exported to other scenarios, changing the way how data is feeded

python 3
"""

import os
import time
import logging
import argparse
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
import fitsio

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s',
                    level=logging.INFO,)

def mask_expnum(arg):
    """ Method to create the mask for each CCD, for the current exposure
    Parameters
    ----------
    exp: int
        Esposure number
    df: dataframe
        Dataframe containing exposure number and the path to catalogue files
        per CCD
    fwhm: str
        Keyword used for FWHM in pixels, column in the catalogue results
    xcentroid: str
        Keyword used for object centroid, in the catalogue results table
    ycentroid: str
        Keyword used for object centroid, in the catalogue results table
    ccd_key: str
        Keyword for CCDNUM in the FITS header
    exp_key: str
        Keyword for EXPNUM in the FITS header
    ccd_dims: tuple
        Dimensions of the CCD (nrows, ncols)
    var: argparse namespace
        Argparse Namespace object. Remember you can create one by hand
    Returns
    -------
    tuple containing exposure number and number of masks for that exposure
    """
    # Uncompress
    exp, df, fwhm, xcentroid, ycentroid, ccd_key, exp_key, ccd_dims, var = arg
    
    logging.info("Working on expnum: {0}".format(exp))
    dfx = df.loc[df["expnum"] == exp]
    # For each of the CCDs: operate and mask
    N_mask = 0
    for idx, row in dfx.iterrows():
        # Open file and get info
        hdu = fitsio.FITS(row["filename"])
        ccdnum = copy.copy(hdu[0].read_header()[ccd_key])
        expnum = copy.copy(hdu[0].read_header()[exp_key])
        if (int(expnum) != exp):
            t_e = "Mismatch in table vs FITS exposure number:"
            t_e += " {0} vs {1}".format(exp, expnum)
            logging.error(t_e)
            exit(1)
        phot = np.copy(hdu[2].read())
        hdu.close()
        
        # Select only objects having a FWHM > 0
        phot = phot[np.where(phot[fwhm] > 0)]

        # Histogram
        if False:
            plt.close("all")
            fig, ax = plt.subplots()
            kw_hist = {"histtype": "stepfilled",
                       "edgecolor": "goldenrod",
                       "color": "lightsteelblue",}
            ax.hist(phot[fwhm], **kw_hist)
            plt.show()
        
        # Look for objects with FWHM higher than some threshold
        if (var.rad is not None):
            sel = phot[np.where(phot[fwhm] >= var.rad)]
        elif (var.p is not None):
            aux_r = np.percentile(phot[fwhm], var.p)
            sel = phot[np.where(phot[fwhm] >= aux_r)] 
        # If there is at least one object obeying the condition, do
        # the masking
        if (sel.size < 1):
            continue
        
        # For the objects obeying the condition, do the masking
        # Locate the centroids and apply a circular mask of n times
        # the FWHM
        min_fwhm = np.min(sel[fwhm])
        msk = np.zeros(ccd_dims).astype(bool)
        for j in range(sel.shape[0]):
            xc, yc = sel[xcentroid][j], sel[ycentroid][j]
            r = var.mf * sel[fwhm][j]
            # Define a open grid, not dense, just 1 row and 1 col
            y_ccd, x_ccd = np.ogrid[:ccd_dims[0], :ccd_dims[1]]
            # Mask
            fcirc = np.power((xc - x_ccd), 2) + np.power((yc - y_ccd), 2)
            tmp_msk = fcirc <= np.power(r, 2)
            # Apply iterative mask to parent mask
            msk = np.logical_or(msk, tmp_msk)
            # Increase exposure number of masks
            N_mask += 1
        # Change to integer, and apply the bit
        msk = msk.astype(int)
        msk[np.where(msk == 1)] = var.bit
        
        # Create a FITS to harbor the data
        outname = "D{0:08}_c{1:02}_maskObj.fits".format(expnum, ccdnum)
        if (var.path is not None):
            outname = os.path.join(var.path, outname)
        fits_m = fitsio.FITS(outname, "rw")
        hlist = [
            {"name": "EXPNUM", "value": expnum, "comment": "",},
            {"name": "CCDNUM", "value": ccdnum, "comment": "",},
            {"name": "COMMENTS", 
             "value": "Masked FWHM >= {0:.3f}".format(min_fwhm), 
             "comment": "",},
        ]
        haux = fitsio.FITSHDR(hlist)
        fits_m.write(msk, header=haux)
        fits_m[-1].write_checksum()
        fits_m[-1].write_history("by Francisco Paz-Chinchon")
        logging.info("Written mask: {0}".format(outname))
    return (expnum, N_mask)

def get_argparse():
    h0 = "Mask by CCD, based on the FWHM of each object and a input threshold"
    h0 += " or statistics. The mask will be applied to objects with"
    h0 += " FWHM>=threshold"
    arg = argparse.ArgumentParser(description=h0)
    h1 = "CSV table containing \'expnum\' and \'filename\' header columns."
    arg.add_argument("csv", help=h1)
    # Excluding group to select the threshold. Either fixed value or 
    # CCD stats
    gexc = arg.add_mutually_exclusive_group()
    h2 = "Threshold. Radius in pix from which to mask objects. Example: 12"
    h2 += " will mask all objects with raius greater or equal than 12 pix" 
    gexc.add_argument("--rad", help=h2, type=float)
    h3 = "Threshold. Percent (range=0, 100) of the total set of radii values,"
    h3 += " from which to mask. Example: 80 will mask all objects with radius"
    h3 += " greater or equal than the value on the 80%% of distribution of"
    h3 += " objects radii."
    gexc.add_argument("--p", help=h3, type=float)
    #
    aux_bit = 4
    h4 = "Bit to be used for masking. Default: {0}".format(aux_bit)
    arg.add_argument("--bit", help=h4, type=int, default=aux_bit)
    h5 = "Directory where to write out the files. Default is current directory"
    arg.add_argument("--path", help=h5)
    h6 = "Filename of 1-column table with subsample of exposures to be used."
    h6 += " File must contain no header. If not given, all exposures will be"
    h6 += " used."
    arg.add_argument("--sub", help=h6)
    aux_mult = 2.
    h7 = "Multiplicative factor to be used to enlarge masking."
    h7 += " Default: {0}".format(aux_mult)
    arg.add_argument("--mf", help=h7, type=float, default=aux_mult)
    #
    res = arg.parse_args()
    return res

def aux_main():    
    var = get_argparse()
    
    # CCD dimensions
    ccd_dims = (4096, 2048)

    # Labels
    fwhm = "FWHM_IMAGE"
    xcentroid = "X_IMAGE"
    ycentroid = "Y_IMAGE"
    ccd_key = "CCDNUM"
    exp_key = "EXPNUM"

    # Open table with files
    df = pd.read_csv(var.csv)
    df.columns = df.columns.map(str.lower)
    
    # Open table with sub-selection of exposures, belonging to the clusters
    # of interest
    if (var.sub is not None):
        sub_arr = np.genfromtxt(var.sub, dtype=int)
        df = df.loc[df["expnum"].isin(sub_arr)]

    # For each exposure, walk through files, no matter the ordering. 
    #
    # Below methods can be all put into an unique method, to use mp
    pool = mp.Pool(processes=mp.cpu_count())
    logging.info("Multiprocessing with {0} threads".format(mp.cpu_count()))
    inarg = [
        (i, df, fwhm, xcentroid, ycentroid, ccd_key, exp_key, ccd_dims, var)
        for i in df["expnum"].unique() 
    ]
    res = pool.map(mask_expnum, inarg)
    pool.close()

    logging.info("Finished, {0}".format(time.ctime()))
    # Save results of number of masks per expnum
    info = pd.DataFrame(res, columns=["expnum", "n_masked"])
    info.to_csv("masked_per_expnum.csv", header=True, index=False)

    return True

if __name__ == "__main__":
    aux_main()
