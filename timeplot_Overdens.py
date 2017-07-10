"""Modules to fit a set of overdense points with ML
"""

import os
import gc
import socket
import pandas as pd
import numpy as np
import datetime
import fitsio
import scipy
import matplotlib.pyplot as plt
from matplotlib.dates import MonthLocator,DayLocator,DateFormatter
from matplotlib.ticker import MultipleLocator
import recarray_tools as rtool
"""testing with sklearn"""
from sklearn.linear_model import Ridge
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
"""
STEPS
1) fit a poly to overdensity
2) distance from every exposure to the fit: cumulative AND median
3) histogram of distances
4) scatter plot between distances and skybrightness/teff/beff
5) how skyfit_qa PCA coeffs relates with the skybrightness/distances?
"""
class Stat():
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


class Fit():
    def __init__(self):
        self.y = np.load("sub_adu.npy")
        self.x = np.load("sub_row.npy")
        #initial and final rows for the issue region 
        self.row_a = 150
        self.row_z = 202

    def poly(self,doPlot=False):
        #restrict to a smaller box
        B1,B2 = self.row_a,self.row_z
        x0 = self.x[np.where(np.logical_and(self.x>=B1*16,self.x<=B2*16))]
        y0 = self.y[np.where(np.logical_and(self.x>=B1*16,self.x<=B2*16))]
        x1 = x0[np.where(np.logical_and(y0>=1.01,y0<=1.042))]
        y1 = y0[np.where(np.logical_and(y0>=1.01,y0<=1.042))]
        #sort
        aux_sort = np.argsort(x1)
        x1 = x1[aux_sort]
        y1 = y1[aux_sort]
        #create matrix version
        x1_aux = x1[:,np.newaxis]
        x1_aux = x1_aux.reshape(-1,1)
        model_fit = []
        for deg in [1,2,3,4,5,6]:
            model = make_pipeline(PolynomialFeatures(deg),Ridge())
            model.fit(x1_aux,y1)
            model_fit.append(model)
            model = None
        if doPlot:
            #to get values: model.predict(x_for_plot)
            plt.scatter(x1,y1,marker=",",color="gray",label="data",s=5)
            colors = ["yellowgreen","magenta","cyan","gold","coral"]
            colors.append("forestgreen")
            x_plot = np.linspace(min(x1),max(x1),1000)
            for ind,m in enumerate(model_fit):
                y_plot = m.predict(x_plot[:,np.newaxis])
                kw = {"color":colors[ind],"lw":1.2*len(colors)-ind}
                kw["label"] = "degree %d"%(ind+1)
                plt.plot(x_plot,y_plot,"-",**kw)
            #all the points
            plt.title("Overdensity fit, for CCD3 boxing")
            plt.legend(loc="lower left")
            kwpdf = {"format":"pdf","dpi":200}
            plt.savefig("poly_fit_{0}.pdf".format(os.getpid()),**kwpdf)
            plt.show()
        sel_fit = model_fit[-1]
        return sel_fit

    def calc_dist(self,f_x,fname_stat,fname_dbinfo):
        """Method to calculate square distance from each CCD3 to the fit. 
        How it relates with surface brightness? and with IMAGE.SKYVAR{A,B}?
        
        statistics array is a structurrred array with cols
        ('med', 'avg', 'med_n', 'avg_n', 'rms_n', 'unc_n', 'mad_n')
        DB info arrays is also structured with cols
        ('nite', 'expnum', 'band', 'ccdnum', 'exptime')
        """
        sumname = "CummSum.npy"
        stt = np.load(fname_stat) 
        db = np.load(fname_dbinfo)
        """improve: concatenate both arrays
        """
        #setup the indices for the overdensity region
        left = np.arange(self.row_a*16,self.row_z*16)[::16]
        idx_left = []
        for L in left:
            #avoid the edge
            idx_left += list(np.arange(L+1,L+8))
        aux_expn = []
        aux_sum = []
        aux_mean = []
        aux_mad = []
        aux_med = []
        aux_rms = []
        aux_kurt = []
        aux_skw = []
        S = Stat()
        #iterate over each CCD
        for ccd in xrange(stt.shape[0]):
            v = []
            for box in idx_left:
                y1 = stt[ccd]["avg_n"][box]
                y2 = f_x.predict(np.array([box])[:,np.newaxis])[0]
                v.append(np.abs(y1-y2))
            v = np.array(v)
            v_aux = v-np.mean(v)
            aux_expn.append(db[ccd]["expnum"][0])
            aux_sum.append(np.sum(v))
            aux_mean.append(np.mean(v))
            aux_med.append(np.median(v))
            aux_rms.append(S.rms(v))
            aux_mad.append(S.mad(v))
            aux_kurt.append(scipy.stats.kurtosis(v_aux,fisher=True,bias=True))
            aux_skw.append(scipy.stats.skew(v_aux))
        #struct array
        tmp = zip(aux_expn,aux_sum,aux_mean,aux_med,aux_kurt,aux_skw,aux_mad,
                aux_rms)
        dt = np.dtype([("expnum","i4"),("sum","f4"),
                    ("avg","f4"),("median","f4"),
                    ("kurt","f4"),("skew","f4"),
                    ("mad","f4"),("rms","f4")])
        X = np.array(tmp,dtype=dt)
        print "Saving results for cumulative sum of distances"
        print sumname
        np.save(sumname,X)
        return True

    def dist_hist(self,npy_fname):
        """Method to create the histogram for the cumulative distances 
        measurements
        """
        data = np.load(npy_fname)
        #sort and drop duplicates
        data = np.sort(data,order="expnum")
        data = rtool.Use.sarr_drop(data,"expnum",posit="first") 
        #
        P = lambda x: np.percentile(data["sum"],x)
        kwh = dict()
        kwh["bins"] = 30
        kwh["histtype"] = "stepfilled"
        kwh["align"] = "mid"
        kwh["log"] = False
        kwh["lw"] = 2
        kwh2 = {"bins":data["sum"].shape[0],"cumulative":True,"normed":False,
            "histtype":"step","color":"red","log":False,"lw":1}
        fig, ax = plt.subplots(1,2,figsize=(8,5))
        #plot
        ax[0].hist(data["sum"],color="royalblue",**kwh)
        ax[0].hist(data["sum"],label="Entire sample",**kwh)
        ax[0].hist(data["sum"],label="Cummulative",**kwh2)
        ax[1].hist(data["sum"],range=(0,P(75)),color="orange",label="75%",
                    **kwh) 
        ax[1].hist(data["sum"],range=(0,P(50)),color="dodgerblue",label="50%",
                    **kwh)
        ax[1].hist(data["sum"],range=(0,P(25)),color="limegreen",label="25%",
                    **kwh)
        #ax[0,1].scatter(data["sum"],data["rms"])
        #ax[1,0].scatter(data["sum"],data["avg"],c="red")
        #ax[1,0].scatter(data["sum"],data["median"],c="blue")
        #ax[1,0].scatter(data["sum"],data["mad"],c="green")
        #ax[1,1].scatter(data["sum"],data["skew"])
        #range, labels
        t = "Histogram of cumulative distances between data and fit,\n"
        t += "for the overdensity region, CCD3 (left amplifier)"
        plt.suptitle(t)
        ax[0].set_xlabel(r"cumulative $D_{<ADU>_{n}}$ (flux)",fontsize=12)
        ax[0].set_ylabel("N",fontsize=12)
        ax[1].set_xlabel(r"cumulative $D_{<ADU>_{n}}$ (flux)",fontsize=12)
        ax[1].set_ylabel("N",fontsize=12)
        ax[0].set_xlim([-100,3000])
        ax[1].set_xlim([5,40])
        #legend
        handles, labels = ax[1].get_legend_handles_labels()
        ax[1].legend(handles,labels,loc="upper right")
        ax[0].legend(loc="center right")
        #spacing
        plt.subplots_adjust(left=0.08,right=0.98)
        plt.savefig("hist_cumSum.pdf",dpi=250,format="pdf")
        plt.show()

    def add_stat(self,npy_fname,skyccd_fname,skyfp_fname):
        """Method to plot additional statistics regarding sky brightness &
        compare with what I already have
        """
        cum = np.load(npy_fname) 
        skyA = np.load(skyfp_fname)
        skyB = np.load(skyccd_fname)
        #sort and remove duplicates
        cum = np.sort(cum,order="expnum")
        skyA = np.sort(skyA,order="expnum")
        skyB = np.sort(skyB,order="expnum")
        cum = rtool.Use.sarr_drop(cum,"expnum",posit="first")
        skyA = rtool.Use.sarr_drop(skyA,"expnum",posit="first")
        skyB = rtool.Use.sarr_drop(skyB,"expnum",posit="first")
        # 
        print "\nFrom cum sum: {0}, {1}".format(cum.dtype.names,cum.shape) 
        print "From sky FP: {0}, {1}".format(skyA.dtype.names,skyA.shape)
        print "From sky per CCD: {0}, {1}".format(skyB.dtype.names,skyB.shape)
        # 
        for i in xrange(cum.shape[0]):
            c1 = cum["expnum"][i] == skyA["expnum"][i]
            c2 = cum["expnum"][i] == skyB["expnum"][i]
            if ~c1 or ~c2:
                logging.error("Arrays doesn't match")
                exit(1)
        #remember left side is AMP-B
        fig,ax = plt.subplots(2,2,figsize=(10,8))
        kw1 = {"s":10,"marker":"."}
        
        if False:
            #using values from the entire focal plane
            ax[0,0].scatter(cum["sum"],skyA["t_eff"],c="royalblue",
                        label="t_eff",**kw1)
            ax[0,0].scatter(cum["sum"],skyA["b_eff"],c="red",
                        label="b_eff",**kw1)
            #
            ax[0,1].scatter(cum["sum"],skyA["skybrightness"],c='forestgreen',
                        label="skybrightness",**kw1)
            #
            ax[1,0].scatter(cum["sum"],skyA["fwhm_asec"],c="goldenrod",
                        label="fwhm [arcsec]",**kw1)
            #
            color = ["royalblue","red","forestgreen","goldenrod"]
            val = ["t_eff","b_eff","skybrightness","fwhm_asec"]
            lab = ["t_eff","b_eff","skybrightness","fwhm [arcsec]"]
            for i in xrange(4):
                ax[1,1].scatter(cum["sum"],skyA[val[i]],c=color[i],
                            label=lab[i],**kw1)
            for subax in [ax[0,0],ax[0,1],ax[1,0],ax[1,1]]:
                subax.set_xlabel(r"cumulative $D_{<ADU>_{n}}$ (flux)")
                subax.legend()
            t = "Cumulative distance vs firstcut_eval parameters\n"
            t += "75% of the sample (3rd quartile)"
            plt.suptitle(t)  
            plt.subplots_adjust(left=0.04,bottom=0.07,right=0.98,top=0.93,
                            wspace=0.21)
            ax[1,1].set_yscale("log")
            ax[0,0].set_xlim([6,37])
            ax[1,0].set_xlim([6,37])
            ax[0,1].set_xlim([6,37])
            ax[1,1].set_xlim([6,37])
            ml1 = MultipleLocator(1)
            for subax in ax.reshape(-1):
                subax.xaxis.set_minor_locator(ml1)
            ax[0,0].yaxis.set_minor_locator(MultipleLocator(0.5))
            ax[0,1].yaxis.set_minor_locator(MultipleLocator(0.1))
            ax[1,0].yaxis.set_minor_locator(MultipleLocator(0.05))
            out_fnm = "cumsum_vs_params1.pdf"
        elif False:
            ax[0,0].scatter(cum["sum"],cum["rms"],c="royalblue",
                            label="rms, 75%",**kw1)
            ax[0,1].scatter(cum["sum"],cum["mad"],c="orange",
                            label="MAD, 75%",**kw1)
            ax[1,0].scatter(cum["sum"],cum["skew"],c="green",
                            label="skewness, 75%",**kw1)
            ax[1,1].scatter(cum["expnum"],cum["sum"],c="goldenrod",**kw1)
            for subax in ax.reshape(-1)[:-1]:
                subax.set_xlabel(r"cumulative $D_{<ADU>_{n}}$ (flux)")
                subax.legend()
            ax[1,1].set_xlabel("expnum")
            ax[1,1].set_ylabel(r"cumulative $D_{<ADU>_{n}}$ (flux)")
            #annotate
            ax[1,1].annotate("75%",xy=(580000,37),xytext=(582500,100),
                            arrowprops=dict(facecolor="lavender",shrink=0.05),)
            ax[0,1].annotate("sinusoidal?",xy=(20,0.01),xytext=(25,0.015),
                            arrowprops=dict(facecolor="darkkhaki",shrink=0.05))
            ax[1,1].axhline(y=37,lw=2,c="blue")
            t = "Cumulative distance statistics"
            plt.suptitle(t)  
            plt.subplots_adjust(left=0.05,bottom=0.13,right=0.98,top=0.93,
                            wspace=0.21)
            #minor ticks 
            for subax in ax.reshape(-1)[:-1]:
                subax.xaxis.set_minor_locator(MultipleLocator(1))
            ax[0,0].yaxis.set_minor_locator(MultipleLocator(0.1))
            ax[0,1].yaxis.set_minor_locator(MultipleLocator(0.001))
            ax[1,0].yaxis.set_minor_locator(MultipleLocator(0.5))
            ax[0,0].set_xlim([6,37])
            ax[1,0].set_xlim([6,37])
            ax[0,1].set_xlim([6,37])
            ax[0,0].set_ylim([0,1.4])
            ax[0,1].set_ylim([0,0.025])
            ax[1,1].set_yscale("log")
            #rotate labels
            plt.setp(ax[1,1].xaxis.get_majorticklabels(),rotation=70)
            out_fnm = "cumsum_vs_params2.pdf"
        else:
            #using values per CCD/AMP
            ax[0,0].scatter(skyB["skyvarb"],cum["sum"],c="mediumturquoise",
                        label="cumulative sum, 75%",**kw1)
            ax[0,1].scatter(skyB["skyvarb"],cum["mad"],c="olive",
                        label="MAD, 75%",**kw1)
            ax[1,0].scatter(skyB["skyvarb"],cum["rms"],c="coral",
                        label="rms, 75%",**kw1)
            ax[1,1].scatter(skyB["expnum"],skyB["skyvarb"],c="royalblue",
                        **kw1)
            for subax in ax.reshape(-1)[:-1]:
                subax.set_xlim([360,590])
                subax.set_xlabel("skyvar ampB")
                subax.xaxis.set_minor_locator(MultipleLocator(10))
                subax.legend()
            ax[1,1].set_xlabel("expnum")
            ax[1,1].set_ylabel("skyvar ampB")
            t = "Sky variance for CCD3, AmpB (issue region) vs"
            t += " cumulative distance parameters\n"
            plt.suptitle(t)
            ax[0,1].annotate("MAD shows a trend",
                            xy=(425,0.0125),xytext=(450,0.02),
                            arrowprops=dict(facecolor="lavender",shrink=0.05))
            #rotate labels
            plt.setp(ax[1,1].xaxis.get_majorticklabels(),rotation=70)
            #minor ticks
            ax[0,0].yaxis.set_minor_locator(MultipleLocator(1))
            ax[0,1].yaxis.set_minor_locator(MultipleLocator(5E-4))
            ax[1,0].yaxis.set_minor_locator(MultipleLocator(0.1))
            ax[1,1].yaxis.set_minor_locator(MultipleLocator(10))
            ax[0,0].set_ylim([6,37])
            ax[0,1].set_ylim([0.003,0.025])
            ax[1,0].set_ylim([0,1.3])
            ax[1,1].set_ylim([360,590])
            plt.subplots_adjust(left=0.04,bottom=0.12,right=0.98,top=0.94,
                            wspace=0.21,hspace=0.24)
            out_fnm = "cumsum_vs_params3.pdf" 

        plt.savefig(out_fnm,dpi=300,format="pdf")
        plt.show()


if __name__ == "__main__":
    print socket.gethostname()
    print "Fit is done using ML, statistics are done for the issue region\n"
    stat_boxes = "stat_ccd3_16x128_uniNmed.npy"
    dbinfo_boxes = "info_ccd3_16x128_uniNmed.npy"
    npy_fname="CummSum.npy"
    sky_ccd = "skybright_ccd.npy"
    sky_fp = "skybright_fp.npy"
    recalculate = False
    plot_hist = False
    add_stat = False
    F = Fit()
    selected_fit = F.poly()
    if recalculate:
        F.calc_dist(selected_fit,stat_boxes,dbinfo_boxes)
    elif plot_hist:
        F.dist_hist(npy_fname)
    elif add_stat:
        F.add_stat(npy_fname,sky_ccd,sky_fp)
    else:
        print "So... no plot?"
        pass
