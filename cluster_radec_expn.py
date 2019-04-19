""" 
Simple clustering in RA-DEC 
"""

import os
import argparse
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Fro gridspec
from matplotlib.gridspec import GridSpec
# For colorbar positioning
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import astropy.units as u
from sklearn.cluster import DBSCAN

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s',
                    level=logging.INFO,)

def sel_by_separation(df, threshold):
    ''' 
    Not optimized at all
    Function to calculate separation assuming plane geometry
    To code it faster, use the separation() method in astropy
    '''
    # List of indices to drop
    idx_drop = []
    for idx1, row1 in df.iterrows():
        for idx2, row2 in df.iterrows():
            if (idx1 == idx2):
                v1 = row1['radec_deg'].separation(row2['radec_deg'])
            else:
                v2 = row1['radec_deg'].separation(row2['radec_deg'])
                # If only one separation is below threshold, then discard 
                # current expnum
                if v2.deg < threshold:
                    idx_drop.append(idx1)
                    break
    try:
        df_sel = df.drop(axis=0, index=idx_drop) 
    except:
        raise
    return df_sel

def diagnostic_plot(skycoo, labels=None, savep=False):
    """ Method to do a quick visualization of the observed fields
    Parameters
    ----------
    skycoo: SkyCoord object
        astropy SkyCoord object contanining the positions to be plotted
    labels: array
        one dimensional array containing the labels for each of the clusters
    savep: boolean
        Wether to save or not the plot
    """
    plt.close("all")
    fig = plt.figure(figsize=(8,6), constrained_layout=True)
    gs = GridSpec(3, 1, figure=fig)
    # Add subplot by the old way
    ax0 = fig.add_subplot(gs[:2, 0], projection="aitoff")
    ax1 = fig.add_subplot(gs[2:, 0])
    # Grid
    ax0.grid(True, which="both", color="gray")
    # RA-DEC
    # We need to use radians, even when the ticks appears in degrees
    aux_ra = skycoo.ra.wrap_at(180 * u.deg).radian
    aux_dec = skycoo.dec.radian
    if (labels is None):
        clr = "blue"
    else:
        clr = labels
    ax0.scatter(aux_ra, aux_dec, marker="h", 
                c=clr, edgecolor="gray", 
                cmap="rainbow", alpha=0.5,)
    # Number of elements per cluster
    ax1.hist(labels, bins=np.unique(labels).size - 1, histtype="stepfilled",
             edgecolor="orange", facecolor="deepskyblue", lw=2)
    ax1.set_xlabel("cluster label")
    ax1.set_ylabel("N")
    #
    plt.subplots_adjust(hspace=0.2, 
                        left=0.07, bottom=0.07, 
                        top=0.99, right=0.99)
    # Save figure
    if savep:
        outnm = "diag_clust_PID{0}.png".format(os.getpid())
        plt.savefig(outnm, dpi=300, format="png")
        logging.info("Written output plot: {0}".format(outnm))
    plt.show()
    return True

def diagnostic_distanceM(distanceMatrix, savep=False):
    """ Simple matrix plotting for quick assessment
    Parameters
    ----------
    distanceMatrix: 2D array
        Square matrix of distances
    savep: boolean
        Whether to save or not the plot
    """
    plt.close("all")
    fig, ax = plt.subplots()
    im = ax.matshow(distanceMatrix)
    # Positioning of the colorbar
    divider = make_axes_locatable(ax)
    caxis = divider.append_axes("right", size="4%", pad=0.3)
    plt.colorbar(im, cax=caxis, extend="max", label=r"dist$_{n}$ [deg]")
    # Easy layout
    plt.tight_layout()
    # Save plot
    if savep:
        outnm = "diag_distMatrix_PID{0}.png".format(os.getpid())
        plt.savefig(outnm, dpi=300, format="png")
        logging.info("Written output plot: {0}".format(outnm))
    plt.show()
    return True

def DBSCAN_radec(skycoo, arg_mdist=None, max_dist=2., min_N=1):
    """ Clustering over RA-DEC coordinates. To do this, distance matrix needs 
    to be calculated.
    Parameters
    ----------
    skycoo: SkyCoord() 
        sky coordinate object from astropy
    arg_mdist: ndarray
        either None to calculate the distance matrix, or a square matrix being
        the distance matrix
    max_dist: float
        maximum distance (in degrees) 2 points should have to be considered 
        as part of the same cluster. The default value=2 is aprox the width 
        of DECam focal plane
    min_N: int
        minimum number of elements a cluster can have
    Return
    ------
    clust: DBSCAN object
        scikit learn DBSCAN object
    """
    # Distance Matrix in sphericals
    # Remember to keep the ordering
    # http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html
    if (arg_mdist is None):
        fd1 = lambda pos_x, pos_all: map(pos_x.separation, pos_all)
        d1 = fd1(skycoo[0], skycoo)
        # Array to harbor the distance matrix
        # With a bit of work this can be parallelized
        logging.warning("Not parallel :|")
        mdist = np.zeros((len(skycoo), len(skycoo)))
        for idx, x in enumerate(skycoo):
            tmp = fd1(x, skycoo)
            tmp = map(lambda y: y.deg, tmp)
            tmp = np.array(list(tmp))
            mdist[idx, :] = tmp
            if not (idx % 10):
                print(idx)
        # Save result 
        try:
            out1 = "mdist.npy"
            np.save(out1, mdist)
            logging.info("Distance matrix saved: {0}".format(out1))
        except:
            raise
    else:
        logging.info("Loading data: {0}".format(arg_mdist))
        mdist = np.load(arg_mdist)
    # Options os clustering where distance is a parameter:
    # - agglomerative clustering
    # - DBSCAN
    # If DBSCAN needs to be called using one of the predefined metrics, use:
    # coo = list(zip(skycoo.ra.deg, skycoo.dec.deg))
    # And then clust.fit(coo)
    clust = DBSCAN(eps=max_dist, 
                   min_samples=min_N, 
                   metric="precomputed")
    clust = clust.fit(mdist)  
    return clust, mdist

def get_argparse():
    h0 = "Simple clustering in RA, DEC. Units for clustering: asec"
    h0 += " NOTE: clustering uses Euclidean distance, be aware of this."
    h0 += " It's possible to improve using spherical coord methods from"
    h0 += "  astropy, introducing it as the metric" 
    arg = argparse.ArgumentParser(description=h0)
    aux_rad = 2.
    h1 = "Diameter in deg to define the clusters. Default: {0}".format(aux_rad)
    arg.add_argument("--rad", help=h1, type=float, default=aux_rad)
    h2 = "Table with ra, dec columns. Units for RA being hh:mm:ss"
    h2 += " and for DEC deg:min:sec"
    arg.add_argument("--tab", help=h2)
    h3 = "Numpy distance matrix, previously calculated. Ordering must be"
    h3 = " the same. Format: npy ndarray"
    arg.add_argument("--mdist", help=h3)
    h4 = "Flag to write plots as PNG files"
    arg.add_argument("--savep", help=h4, action="store_true")
    arg = arg.parse_args()
    return arg

def aux_main():
    # Get the arguments
    var = get_argparse()
    # Read table with header containing ra, dec (can contain more)
    df = pd.read_csv(var.tab)
    df.columns = df.columns.map(str.lower)
    if (("ra" not in df.columns) or ("dec" not in df.columns)):
        logging.error("Column headers missing: ra, dec")
        exit(1)
    # FK5 is the coordinate system for the telescope
    logging.info("Using FK5 as the telescope coordinate system")
    fk5c = SkyCoord(df['ra'], df['dec'], 
                    unit=(u.hourangle, u.deg),
                    frame=FK5)
    # DBSCAN
    clust, mdist = DBSCAN_radec(fk5c, arg_mdist=var.mdist, 
                                max_dist=var.rad, min_N=1)
    
    # Add new column to the input table, where the cluster ID is written
    df["cluster_id"] = clust.labels_
    outc = "cluster_res_pid{0}.csv".format(os.getpid())
    df.to_csv(outc, index=False, header=True)
    logging.info("Saved: {0}".format(outc))
    
    # Diagnostic of exposures in the Sky
    diagnostic_plot(fk5c, labels=clust.labels_, savep=var.savep)
    # Diagnostic plot distance matrix
    diagnostic_distanceM(mdist, savep=var.savep)
    return True

if __name__ == '__main__':

    aux_main()
