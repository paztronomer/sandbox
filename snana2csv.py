
# coding: utf-8

# In[1]:

import pandas as pd
import numpy as np
import os 
import socket
import vaex
import matplotlib.pyplot as plt


# In[40]:

# Construct a single DF with all the lightcurves
def snana_fmt(fname):
    x = []
    with open(fname, "r") as g:
        for line in g:
            if (line[:5] == "SNID:"):
                id = line[5:].strip()
            if (line[:4] == "OBS:"):
                aux = line[4:].split()
                x.append([id] + aux)
    return x

def df_constr(data):
    # Receives a list os lists
    cols = ["SNID", "MJD","FLT","FIELD","FLUXCAL","FLUXCALERR","PHOTFLAG","PHOTPROB"]
    cols += ["ZPFLUX","PSF","SKYSIG","SKYSIG_T","GAIN","XPIX","YPIX","NITE"]
    cols += ["EXPNUM","CCDNUM","OBJID"]
    d = dict(zip(cols, data))
    # Change data type
    try:
        for key in ["SNID", "NITE", "EXPNUM", "CCDNUM", "OBJID"]:
            d[key] = map(int, d[key])
        for key_float in ["MJD", "PHOTFLAG", "FLUXCAL", "FLUXCALERR", "PHOTPROB", "ZPFLUX", 
                          "PSF", "SKYSIG", "SKYSIG_T", "GAIN", "XPIX", "YPIX"]:
            d[key_float] = map(float, d[key_float])
        # Re-make the PHOTFLAG
        d["PHOTFLAG"] = map(int, d["PHOTFLAG"])
    except:
        print data
    df = pd.DataFrame(d)
    print df.info()
    # Write out
    outfnm = os.path.join(os.getcwd(), "LC_output_IC170922_dp801_20171007.csv")
    df.to_csv(outfnm, index=False, header=True)
    return True
    
def walk(pth):
    store = []
    for root, dirs, files in os.walk(pth):
        for f in files:
            if ".dat" in f:
                tab = snana_fmt(os.path.join(root, f))
                store += tab
    # Create the DF
    tmp_store = map(list, zip(*store))
    df_constr(tmp_store)
    
def quick_show(fnm1, fnm2): 
    # Load dataframe
    df1 = pd.read_csv(fnm1)
    df2 = pd.read_csv(fnm2)
    # ds = vaex.from_pandas(df, name="lc")
    plt.close("all")
    plt.scatter(df1["MJD"], (df1["FLUXCAL"]), c=df1["SKYSIG"], marker=".", 
                cmap="gist_rainbow")
    plt.scatter(df2["MJD"], (df2["FLUXCAL"]), c=df2["SKYSIG"], marker=".", 
                cmap="gist_rainbow")
    # ds.plot("MJD", "FLUXCAL")
    print df1.info(), df2.info()
    plt.show()


# In[41]:

if __name__ == "__main__":
    print socket.gethostname()
    if False:
        # To construct the DF from SNANA files
        # output_IC170922_dp801_20171007
        # output_IC170922_dp802_20171009
        folder = os.path.join(os.getcwd(), "output_IC170922_dp801_20171007")
        walk(folder)
    # Make the plots
    csv1 = "LC_output_IC170922_dp801_20171007.csv"
    csv2 = "LC_output_IC170922_dp802_20171009.csv"
    quick_show(csv1, csv2)


# In[ ]:




# In[ ]:



