""" Simple script to impose a threshold in the RA-DEC ponting of the telescope
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import astropy.units as u

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

if __name__ == '__main__':
    # Load table
    fnm = 'r4133_nobjects.csv'
    fnm_out = 'r4133_sel_by_separation.csv'
    
    print('Loading: {0}'.format(fnm))
    
    df = pd.read_csv(fnm)
    df.columns = df.columns.map(str.lower)
    # Transform TELRA/TELDEC to ra/dec in deg 
    fk5c = SkyCoord(df['telra'], df['teldec'], unit=(u.hourangle, u.deg),
                    frame=FK5)
    # Pass this coordinate object to pandas df
    df['radec_deg'] = fk5c
    
    minsep = 30 / 3600.# asec in deg
    
    # Apply a cut learn from the plot of nobjects vs exptime
    df = df.loc[df['nobjects'] < 100000]
    df_sel = sel_by_separation(df, minsep)

    # Write out
    print('Saving results: {0}'.format(fnm_out))
    df_sel.to_csv(fnm_out, index=False, header=True)
