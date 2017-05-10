import os
import time
import numpy as np
import tables
import matplotlib
import matplotlib.pyplot as plt

class View():
    def level_view(self,level,label=None):
        """Receives a list containingthe arrays for each decomposition:
        approximate, horizontal, vertical, and diagonal. Plots the level.
        Inputs:
        - level: list of 2D arrays, on which each one is a axial decomposition,
        being cA, cH, cV, cD
        """
        #note "axy" is a ndarray on which each element is a axes object.
        #for sharex/y other values are row, col
        fig,axy = plt.subplots(2,2,figsize=(7,7),sharex="all",sharey="all")
        norm = matplotlib.colors.Normalize()
        kw = {"origin":"lower","interpolation":"none","cmap":"viridis",
            "norm":norm}
        print len(level)
        for x in level:
            print x.shape
            print x.min()
            print x.max()
        extr = [np.array(np.min(x),np.max(x)) for x in level]
        print extr
        extr = [np.min(extr),np.max(extr)]
        print extr
        axy[0,0].imshow(level[0],**kw)
        axy[0,1].imshow(level[1],**kw)
        axy[1,0].imshow(level[2],**kw)
        axy[1,1].imshow(level[3],**kw)
        #hide ticks from right column plots
        plt.setp([a.get_xticklabels() for a in axy[0,:]],visible=False)
        plt.setp([a.get_yticklabels() for a in axy[:,1]],visible=False)
        if ~(label is None): fig.suptitle(label)
        fig.subplots_adjust(hspace=0.01,wspace=0.01)
        plt.show()
        fig,axy = None,None
        return True

    def levels(self,h5tab,fname=None):
        """Method to plot different levels of the DWT
        """
        col = h5tab.colnames
        Nrow = h5tab.nrows/len(col)
        binned = fname[fname.find("_b")+1:fname.find(".h5")]
        #plot every decomposition level
        for cc in col:
            #L = [rr[cc] for rr in h5tab.iterrows()]
            L = [rr[cc] for rr in h5tab.iterrows()]
            print cc,len(L)
            auxlabel = "Level {0}, {1} ({2})".format(cc,binned,h5tab.name)
            View().level_view(L,label=auxlabel)
        #http://www.pytables.org/usersguide/libref/structured_storage.html?
        #highlight=colnames#tables.Table.colnames
        print h5tab.name
        print h5tab.title
        print h5tab.colnames
        print 'Number os columns: ',h5tab.cols
        print 'Number of rows: ',h5tab.nrows
        print 'Shape: ',h5tab.shape
        print len(h5tab.iterrows()[:])
        #for i in xrange()
        #[row["c_A"] for row in h5table.iterrows()]


class OpenH5():
    """class oriented to open and close a pytables instance. It"s not the same
    to close the instance here inside the same class than outside
    """
    def __init__(self,folder,fname):
        fname = os.path.join(folder,fname)
        self.h5 = tables.open_file(fname,driver="H5FD_CORE")

    def closetab(self):
        self.h5.close()


if __name__ == "__main__":
    wavelet = "dmey"
    folder = "/Users/fco/Code/des_calibrations/dwt_files"
    folder = os.path.join(folder,wavelet)
    #to retrict the depth to be walked!
    DEPTH = 0
    for root,dirs,files in os.walk(folder):
        if root.count(os.sep) >= DEPTH:
            del dirs[:]
        #folder = "/Users/fco/Code/des_calibrations/dwt_files/haar"
        #haar_feat1_b6464.h5
        FPname = [binned for binned in files if (".h5" in binned)]
        for aux in FPname:
            print "\n\t{0}".format(aux)
            h5aux = OpenH5(folder,aux).h5
            if wavelet == "dmey":
                h5tab = h5aux.root.dflat.dmey
            elif wavelet == "haar":
                h5tab = h5aux.root.dflat.haar
            View().levels(h5tab,fname=aux)
