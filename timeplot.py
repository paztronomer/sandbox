"""Modifiable script to time series plotting
"""
import sys
import os
import gc
import socket
import pandas as pd
import numpy as np
from numpy import ma
import datetime
import fitsio
import matplotlib.pyplot as plt
from matplotlib.dates import MonthLocator, DayLocator, DateFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from sklearn.linear_model import Ridge
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline


class FPBinned():
    def __init__(self, fpath):
        f = fitsio.FITS(fpath)
        # x_header = fitsio.read_header(fpath)
        x_head = f[0].read_header()
        x_hdu = f[0].read()
        self.ccd = np.copy(x_hdu)
        self.header = x_head
        f.close()
        # do I need to close the fits? M_hdu.close()


class Tide():
    def __init__(self, stat=None, db=None):
        self.stt = np.load(stat)
        self.db = np.load(db)

    def ccd_x(self, pre_nm="plt_16x128_NperCCD"):
        """Ad-hoc function to plot a subset of the array data
        1) The structure of the statistics array (a structurred one) is: 2
        dimensions, (N1, N2) where N2 is the amount of subsections
        (increasing in dim0 and inside in dim1) and N1 is the number
        of analyzed CCDs. For each one of the N2 entries, the columns are:
        ('med', 'avg', 'med_n', 'avg_n', 'rms_n', 'unc_n', 'mad_n')
        2) for the database info array, the dimensions are (N1, 1) where N1 is
        the number of analyzed CCDs. The columns are:
        ('nite', 'expnum', 'band', 'ccdnum', 'exptime')
        """
        print self.stt.dtype.names
        # First:
        # setup the indices belonging to the issue region
        # indices for the region of the anomaly dim0:2208-3248
        ccd_l = np.arange(2208, 3248)[::16]
        idx_left = []
        for x in ccd_l:
            idx_left += list(np.arange(x, x+9))
        idx_left = np.array(idx_left)
        # indices for the region at the edges
        edge_l = np.arange(0, 4096)[::16]
        idx_edge = []
        for idx, x in enumerate(edge_l):
            if (idx > 0): 
                idx_edge += [x-1]
            idx_edge += [x]
        idx_edge = np.array(idx_edge)
        # With the followig, I'll remove the border regions. This works
        # for 1D arrays
        aux0 = np.arange(0, idx_left.shape[0])
        # returns the indices where the elements of the second array must
        # be inserted to preserve the order
        aux1 = np.searchsorted(idx_left, idx_edge, side="left")
        # returns a new array from the first, discarding the elements whose
        # indices are in the second array
        aux2 = np.delete(aux0, aux1)
        idx_left = idx_left[aux2]
        # indices not in the left side neither on the edges
        Ndim0, Ndim1 = 256, 16
        oth0 = np.arange(0, Ndim0*Ndim1)
        oth1 = oth0[np.delete(
                np.arange(0, oth0.shape[0]), 
                np.searchsorted(oth0, np.array(idx_left))
                )]
        idx_other = oth1[np.delete(
                        np.arange(0, oth1.shape[0]), 
                        np.searchsorted(oth1, np.array(idx_edge))
                        )]
        # iteratively fill the plot
        #
        # 1) define the locators for date-axis
        months = MonthLocator()
        days = DayLocator()
        dateFmt = DateFormatter("%b-%d-%Y")# b if for Month name
        months2 = MonthLocator()
        days2 = DayLocator()
        dateFmt2 = DateFormatter("%b-%d-%Y")
        # transform nite to a plottable quantity
        xdate = [datetime.datetime.strptime(str(date[0]), "%Y%m%d")
                for date in self.db[:]["nite"]]
        #
        # 2) setup plot
        gc.collect()
        plt.close("all")
        fig = plt.figure(figsize=(9, 8))
        ax1 = fig.add_subplot(211, axisbg="1.")
        ax2 = fig.add_subplot(212, axisbg="1.")
        karg = {"marker":"o", "s":20, "alpha":1, "edgecolor":"black"}
        karg0 = {"marker":".", "s":10, "alpha":1, "c":"red"}
        karg0.update({"edgecolor":"none"})
        karg1 = {"marker":".", "s":10, "alpha":1, "c":"black"} 
        karg1.update({"edgecolor":"none"})
        karg2 = {"marker":".", "s":100, "alpha":.3, "c":"forestgreen"}
        karg3 = {"marker":".", "s":100, "alpha":.3, "c":"blue"}
        #
        # 3) putt the points where they belongs
        # iterate over different expnum
        for i in np.arange(self.stt.shape[0]):# [::6]:
            print i
            gc.collect()
            """use avg_n, rms_n"""
            key1 = "avg_n"
            key2 = "rms_n"
            """plot issue region"""
            x0 = np.array([xdate[i]]*idx_left.shape[0])
            y1 = np.take(self.stt[i][key1], idx_left)
            y2 = np.take(self.stt[i][key2], idx_left)
            ax1.scatter(x0, y1, **karg2)
            ax2.scatter(x0, y2, **karg3)
            """plot all the other boxes"""
            aux_date = np.array([xdate[i]]*self.stt.shape[1])
            im1 = ax1.scatter(aux_date, self.stt[i][key1], **karg1)
            im2 = ax2.scatter(aux_date, self.stt[i][key2], **karg1)
            """plot the vertical edges of the CCD"""
            ed_x0 = np.array([xdate[i]]*idx_edge.shape[0])
            ed_y1 = np.take(self.stt[i][key1], idx_edge)
            ed_y2 = np.take(self.stt[i][key2], idx_edge)
            ax1.scatter(ed_x0, ed_y1, **karg0)
            ax2.scatter(ed_x0, ed_y2, **karg0)
        #
        # 4) post setup
        ax1.set_ylabel(r"$\log_{10} \langle x \rangle_{norm}$", fontsize=14)
        ax2.set_ylabel(r"$\log_{10}$ RMS$_{norm}$", fontsize=14)
        # cb1 = fig.colorbar(im1, ax=ax1)
        # cb2 = fig.colorbar(im2, ax=ax2)
        # cb1.set_label(r"$\langle x \rangle$")
        # cb2.set_label(r"$\langle x \rangle$")
        # ax1.set_ylim([-100, 400])
        ax1.xaxis.set_major_locator(months)
        ax1.xaxis.set_major_formatter(dateFmt)
        ax1.xaxis.set_minor_locator(days)
        ax1.autoscale_view()
        ax1.grid(True)
        ax2.xaxis.set_major_locator(months2)
        ax2.xaxis.set_major_formatter(dateFmt2)
        ax2.xaxis.set_minor_locator(days2)
        ax2.autoscale_view()
        ax2.grid(True)
        fig.autofmt_xdate()
        ax2.set_yscale("log")
        ax1.set_yscale("log")
        ax1.autoscale_view()
        ax2.autoscale_view()
        plt.subplots_adjust(left=0.1, bottom=0.1, top=0.98, right=0.99, 
                        hspace=0.10)
        if True:
            outnm = "{0}_pid{1}.png".format(pre_nm, os.getpid())
            plt.savefig(outnm, dpi=400, facecolor="w", 
                    edgecolor="w", 
                    orientation="portrait", papertype=None, format="png", 
                    transparent=False, bbox_inches=None, pad_inches=0.1, 
                    frameon=None)
        gc.collect()
        plt.show()

    def row_profile(self, ccd, amp=0, suffix=None, saveplot=True, 
                NCCD=623, NBOX=4096):
        """ Method for plotting the boxes subdividing the CCD, focused
        on one amplifier only (the one containing the issue shown in the
        skytemplates)
        Inputs
        - ccd: integer describing the ccd (0-99)
        - amp: integer describing the amplifier to be used. If 0, use the left
        amplifier, 1 for the right amplifier
        - suffix: string to be added to the output filename
        - saveplot: wheter to save or not the plot
        - NCCD: number of CCDs present in the timeseries array
        - NBOX: number of boxes subdividing each of the CCDs
        """
        print "Statistics array columns: {0}".format(self.stt.dtype.names)
        print "DB array columns: {0}".format(self.db.dtype.names)
        # boxing of 16x128 produces 256 rows and 16 columns
        row_ini, row_end = 0, 256
        row_issue1, row_issue2 = 150, 202
        if amp == 0:
            # left amplifier without the edge box
            left = np.arange(row_ini * 16, row_end * 16)[::16]
            idx_left = []
            for L in left:
                # avoid the edge
                idx_left += list(np.arange(L + 1, L + 8))
            idx = np.array(idx_left)
        elif amp == 1:
            # right amplifier without the edge box
            right = np.arange(row_ini*16 + 8, row_end*16 + 8)[::16]
            idx_right = []
            for R in right:
                # avoid the edge
                idx_right += list(np.arange(R, R + 8))
            idx = np.array(idx_right)
        else:
            logging.error("Error in amplifier selection")
            exit(1)

        plt.close("all")
        fig, ax = plt.subplots(2, 2, figsize=(13,9)) 
        # Auxiliary masked array for plotting row numbes
        row_mask_aux = np.ones((NCCD, NBOX), dtype=bool)
        row_mask_aux[:, idx] = False
        row_plot = np.tile(np.arange(NBOX), NCCD).reshape(NCCD, NBOX)
        row_plot = ma.array(row_plot, mask=row_mask_aux)
        # For the boxes, I need to select only left side, using a mask
        # Apply the mask to the boxes-statistics and db-info arrays
        not_idx = np.ones((NCCD, NBOX),dtype=bool)
        not_idx[:, idx] = False
        stt_sub = ma.array(self.stt, mask=not_idx)
        #
        # 1st subplot, simple profile
        # Use colormap for the value used on normalization
        kw_ax0 = {
            "c": stt_sub["norm"].compressed().ravel(),
            "s": 10, 
            "marker": ".", 
            "alpha": 0.8, 
            "edgecolors":"none",
            "cmap": "plasma",}
        ax[0, 0].axvspan(row_issue1 * 16, row_issue2 * 16, 
                        facecolor="lightgray", alpha=0.5)
        p0 = ax[0, 0].scatter(row_plot.compressed().ravel(), 
                            stt_sub["avg_n"].compressed().ravel(), 
                            **kw_ax0)
        ax[0, 0].set_xlabel("Row (box) number")
        ax[0, 0].set_ylabel(r"Boxes values, $\langle ADU \rangle_{norm}$")
        ax[0, 0].set_title("Row profile")
        cb0 = plt.colorbar(p0, ax=ax[0, 0])
        cb0.set_label("Normalization value (median)")
        majorLocator_x0 = MultipleLocator(512)
        majorFormatter_x0 = FormatStrFormatter("%d")
        minorLocator_x0 = MultipleLocator(128)
        ax[0, 0].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 0].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 0].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 0].set_xlim([0,NBOX])
        # For CCD2, AmpA ylim [-0.00118, -0.0006]
        # For CCD3, AmpB ylim [-30,220]
        ax[0, 0].set_ylim([-30,220]) 
        #
        # 2nd value of boxes vs nite, with quartiles overplotted for each nite
        # Define the date formatter for the x axis
        months = MonthLocator()
        days = DayLocator()
        dateFmt = DateFormatter("%b/%d/%y")
        months2 = MonthLocator()
        days2 = DayLocator()
        dateFmt2 = DateFormatter("%b/%d/%y")
        # Sort by expnum and transform nite to a plottable quantity
        xdate = [datetime.datetime.strptime(str(date[0]), "%Y%m%d")
                for date in np.sort(self.db["nite"])]
        kw_ax1 = {
            "s": 1,
            "c": "silver",
            "marker": ".",
            }
        q0, q1, q2, q3, q4 = [], [], [], [], []
        # For each nite, plot profile
        for n in xrange(len(xdate)):
            # x.compressed() is equivalent to x[~x.mask]
            tmp_stt = stt_sub["avg_n"][n, :].compressed()
            tmp_date = [xdate[n]] * tmp_stt.shape[0]
            q0.append(np.percentile(tmp_stt, 5)) 
            q1.append(np.percentile(tmp_stt, 25)) 
            q2.append(np.percentile(tmp_stt, 50))
            q3.append(np.percentile(tmp_stt, 75))
            q4.append(np.percentile(tmp_stt, 95)) 
            ax[0, 1].scatter(tmp_date, tmp_stt, **kw_ax1)
        ax[0, 1].scatter(xdate, q0, c="violet", s=2, marker="^", label="5%")
        ax[0, 1].scatter(xdate, q1, c="red", s=2, marker="s", label="25%")
        ax[0, 1].scatter(xdate, q2, c="blue", s=2, marker="D", label="50%")
        ax[0, 1].scatter(xdate, q3, c="green", s=2, marker="o", label="75%")
        ax[0, 1].scatter(xdate, q4, c="black", s=2, marker="v", label="95%")
        ax[0, 1].set_ylabel(r"Box values, $\langle ADU \rangle_{norm}$")
        ax[0, 1].set_title("Amplifier boxes values vs nite")
        ax[0, 1].xaxis.set_major_locator(months)
        ax[0, 1].xaxis.set_major_formatter(dateFmt)
        ax[0, 1].xaxis.set_minor_locator(days)
        ax[0, 1].autoscale_view()
        kw_grid1 = {
            "color": "lightgray",
            "linestyle": "dotted",
            "dash_capstyle": "round",
            "alpha": 0.7,}
        ax[0, 1].grid(**kw_grid1)
        xlabels_ax1 = ax[0, 1].get_xticklabels()
        plt.setp(xlabels_ax1, rotation=30, fontsize=10)
        # Setup legend
        handles_ax1, labels_ax1 = ax[0, 1].get_legend_handles_labels()
        kw_lab1 = {
            "loc": "lower center",
            "ncol": len(labels_ax1),
            "markerscale": 2.1,
            "frameon": True,
            "fancybox": True,
            "framealpha": 0.5,
            "fontsize": 9,}
        ax[0, 1].legend(handles_ax1, labels_ax1, **kw_lab1)
        # For CCD2, AmpA, ylim [-0.00118, -0.0006]
        # For CCD3, AmpB, ylim [0.88, 1.3]
        ax[0, 1].set_ylim([0.88, 1.3])
        #
        # 3rd density map (small scale when needed)
        kw_ax2 = {
            "bins": 80,
            "cmap": "viridis",
            "cmin": 1,}
        # For CCD3, AmpB must restrict histogram borders
        kw_ax2.update({"range": [[0, NBOX], [-30, 220.]], "bins": 100,})
        ax[1, 0].axvspan(row_issue1 * 16, row_issue2 * 16, 
                        facecolor="lightgray", alpha=0.5)
        p2 = ax[1, 0].hist2d(row_plot.compressed().ravel(), 
                            stt_sub["avg_n"].compressed().ravel(),
                            **kw_ax2)  
        ax[1, 0].set_xlabel("Row (box) number")
        ax[1, 0].set_ylabel(r"Boxes values, $\langle ADU \rangle_{norm}$")
        ax[1, 0].set_title("Density map for row profile".format(ccd))
        cb2 = plt.colorbar(p2[-1], ax=ax[1, 0])
        cb2.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb2.set_label("N")
        majorLocator_x2 = MultipleLocator(512)
        majorFormatter_x2 = FormatStrFormatter("%d")
        minorLocator_x2 = MultipleLocator(128)
        ax[1, 0].xaxis.set_minor_locator(minorLocator_x2)
        ax[1, 0].xaxis.set_major_locator(majorLocator_x2)
        ax[1, 0].xaxis.set_major_formatter(majorFormatter_x2)
        ax[1, 0].set_xlim([0,NBOX])
        #
        # 4th boxes for the issue region versus nite
        months = MonthLocator()
        days = DayLocator()
        dateFmt = DateFormatter("%b/%d/%y")
        months2 = MonthLocator()
        days2 = DayLocator()
        dateFmt2 = DateFormatter("%b/%d/%y")
        # Sort by expnum and transform nite to a plottable quantity
        xdate3 = [datetime.datetime.strptime(str(date[0]), "%Y%m%d")
                for date in np.sort(self.db["nite"])]
        kw_ax3 = {
            "s": 1,
            "c": "lightsteelblue",
            "marker": ".",
            }
        q0, q1, q2, q3, q4 = [], [], [], [], []
        # For each nite, plot profile
        for n in xrange(len(xdate3)):
            # x.compressed() is equivalent to x[~x.mask]
            tmp_stt = stt_sub["avg_n"][n,:].compressed()
            # To select only the issue region, must consider that one box
            # has been removed from the edge, and thus only 7 boxes remains
            tmp_stt = tmp_stt[row_issue1 * 7: row_issue2*7 + 1]
            tmp_date = [xdate3[n]] * tmp_stt.shape[0]
            q0.append(np.percentile(tmp_stt, 5)) 
            q1.append(np.percentile(tmp_stt, 25)) 
            q2.append(np.percentile(tmp_stt, 50))
            q3.append(np.percentile(tmp_stt, 75))
            q4.append(np.percentile(tmp_stt, 95)) 
            ax[1, 1].scatter(tmp_date, tmp_stt, **kw_ax3)
        ax[1, 1].scatter(xdate3, q0, c="violet", s=2, marker="^", label="5%")
        ax[1, 1].scatter(xdate3, q1, c="red", s=2, marker="s", label="25%")
        ax[1, 1].scatter(xdate3, q2, c="blue", s=2, marker="D", label="50%")
        ax[1, 1].scatter(xdate3, q3, c="green", s=2, marker="o", label="75%")
        ax[1, 1].scatter(xdate3, q4, c="black", s=2, marker="v", label="95%")
        ax[1, 1].set_ylabel(r"Box values, $\langle ADU \rangle_{norm}$")
        ax[1, 1].set_title("Region of the issue (ampB) versus nite")
        ax[1, 1].xaxis.set_major_locator(months)
        ax[1, 1].xaxis.set_major_formatter(dateFmt)
        ax[1, 1].xaxis.set_minor_locator(days)
        ax[1, 1].autoscale_view()
        kw_grid3 = {
            "color": "lightgray",
            "linestyle": "dotted",
            "dash_capstyle": "round",
            "alpha": 0.7,}
        ax[1, 1].grid(**kw_grid3)
        # For CCD2, AmpA, ylim [-0.001, -0.00085]
        # For CCD3, AmpB, ylim [0.88, 1.3]
        ax[1, 1].set_ylim([0.88, 1.3])
        # Setup xtick labels 
        xlabels_ax3 = ax[1, 1].get_xticklabels()
        plt.setp(xlabels_ax3, rotation=30, fontsize=10)
        # Setup legend
        handles_ax3, labels_ax3 = ax[1, 1].get_legend_handles_labels()
        kw_lab3 = {
            "loc": "lower center",
            "ncol": len(labels_ax3),
            "markerscale": 2.1,
            "frameon": True,
            "fancybox": True,
            "framealpha": 0.5,
            "fontsize": 9,}
        ax[1, 1].legend(handles_ax3, labels_ax3, **kw_lab3)
        # Final 
        supt = "CCD{0}, amplifier showing the 'specter'.".format(ccd)
        supt += " Boxes of h=16 pix by w=128 pix"
        supt += "\nNOTE: Reduced images (used for skytemplate)."
        plt.suptitle(supt, color="blue")
        plt.subplots_adjust(left=0.1, bottom=0.08, top=0.9, right=0.98, 
                        hspace=0.26, wspace=0.26)
        if suffix is None:
            suffix = "pid" + str(os.getpid())
        if saveplot:
            outname = "rowProfile_{0:02}_{1}.png".format(ccd,suffix)
            plt.savefig(outname, dpi=400, facecolor="w", edgecolor="w", 
                    orientation="portrait", papertype=None, format="png", 
                    transparent=False, bbox_inches=None, pad_inches=0.1, 
                    frameon=None)
        plt.show()

    def zoom_row_profile(self, ccd, amp=0, suffix=None, saveplot=True, 
                NCCD=623, NBOX=4096):
        """ Method based on row_profile(). The difference is this one
        makes a zoom of the interesting zones, discarding ranges of values
        outside this region.
        It plots the boxes subdividing the CCD, focused on the amplifier 
        containig the issue (shown in the skytemplates)
        Inputs
        - ccd: integer describing the ccd (0-99)
        - amp: integer describing the amplifier to be used. If 0, use the left
        amplifier, 1 for the right amplifier
        - suffix: string to be added to the output filename
        - saveplot: wheter to save or not the plot
        - NCCD: number of CCDs present in the timeseries array
        - NBOX: number of boxes subdividing each of the CCDs
        """
        print "Statistics array columns: {0}".format(self.stt.dtype.names)
        print "DB array columns: {0}".format(self.db.dtype.names)
        # boxing of 16x128 produces 256 rows and 16 columns
        row_ini, row_end = 0, 256
        row_issue1, row_issue2 = 150, 202
        if amp == 0:
            # left amplifier without the edge box
            left = np.arange(row_ini * 16, row_end * 16)[::16]
            idx_left = []
            for L in left:
                # avoid the edge
                idx_left += list(np.arange(L + 1, L + 8))
            idx = np.array(idx_left)
        elif amp == 1:
            # right amplifier without the edge box
            right = np.arange(row_ini*16 + 8, row_end*16 + 8)[::16]
            idx_right = []
            for R in right:
                # avoid the edge
                idx_right += list(np.arange(R, R + 8))
            idx = np.array(idx_right)
        else:
            logging.error("Error in amplifier selection")
            exit(1)
        #
        plt.close("all")
        fig, ax = plt.subplots(2, 2, figsize=(13, 9)) 
        # Auxiliary masked array for plotting row numbes
        row_mask_aux = np.ones((NCCD, NBOX), dtype=bool)
        row_mask_aux[:, idx] = False
        row_plot = np.tile(np.arange(NBOX), NCCD).reshape(NCCD, NBOX)
        row_plot = ma.array(row_plot, mask=row_mask_aux)
        # For the boxes, I need to select only left side, using a mask
        # Apply the mask to the boxes-statistics and db-info arrays
        not_idx = np.ones((NCCD, NBOX),dtype=bool)
        not_idx[:, idx] = False
        stt_sub = ma.array(self.stt, mask=not_idx)
        #
        # 1st plot 2D histogram of the entire set of rows
        kw_ax0 = {
            "bins": 100,
            "cmap": "viridis",
            "cmin": 1,}
        # For CCD3, AmpB, dict needs to be updated
        # using ({"range": [[0, NBOX], [-300, 400]],})
        # kw_ax0.update({"range": [[0, NBOX], [-200, 250]],})
        ax[0, 0].axvspan(row_issue1 * 16, row_issue2 * 16, 
                        facecolor="lightgray", alpha=0.5)
        p0 = ax[0, 0].hist2d(row_plot.compressed().ravel(), 
                            stt_sub["avg_n"].compressed().ravel(),
                            **kw_ax0)  
        ax[0, 0].set_xlabel("Row (box) number")
        ax[0, 0].set_ylabel(r"Boxes values, $\langle ADU \rangle_{norm}$")
        ax[0, 0].set_title("Density map for row profile".format(ccd))
        cb0 = plt.colorbar(p0[-1], ax=ax[0, 0])
        cb0.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb0.set_label("N")
        majorLocator_x0 = MultipleLocator(512)
        majorFormatter_x0 = FormatStrFormatter("%d")
        minorLocator_x0 = MultipleLocator(128)
        ax[0, 0].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 0].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 0].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 0].set_xlim([0,NBOX])
        # For CCD2, AmpA, ylim [-0.00118, -0.0006]
        # For CCD3, AmpA, ylim [-25, 160]
        # For CCD2, AmpA, ylim [-140, 180]
        # For CCD3, AmpB, ylim [-17, 140]
        # For CCD2, AmpB, ylim [-215000, 20000]
        ax[0, 0].set_ylim([-215000, 20000])
        #
        # 2nd value of boxes vs nite, with quartiles overplotted for each nite
        # Define the date formatter for the x axis
        months = MonthLocator()
        days = DayLocator()
        dateFmt = DateFormatter("%b/%d/%y")
        months2 = MonthLocator()
        days2 = DayLocator()
        dateFmt2 = DateFormatter("%b/%d/%y")
        # Sort by expnum and transform nite to a plottable quantity
        xdate = [datetime.datetime.strptime(str(date[0]), "%Y%m%d")
                for date in np.sort(self.db["nite"])]
        kw_ax1 = {
            "s": 1,
            "c": "silver",
            "marker": ".",
            }
        q0, q1, q2, q3, q4 = [], [], [], [], []
        # For each nite, plot profile
        for n in xrange(len(xdate)):
            # x.compressed() is equivalent to x[~x.mask]
            tmp_stt = stt_sub["avg_n"][n, :].compressed()
            tmp_date = [xdate[n]] * tmp_stt.shape[0]
            q0.append(np.percentile(tmp_stt, 5)) 
            q1.append(np.percentile(tmp_stt, 25)) 
            q2.append(np.percentile(tmp_stt, 50))
            q3.append(np.percentile(tmp_stt, 75))
            q4.append(np.percentile(tmp_stt, 95)) 
            ax[0, 1].scatter(tmp_date, tmp_stt, **kw_ax1)
        ax[0, 1].scatter(xdate, q0, c="violet", s=2, marker="^", label="5%")
        ax[0, 1].scatter(xdate, q1, c="red", s=2, marker="s", label="25%")
        ax[0, 1].scatter(xdate, q2, c="blue", s=2, marker="D", label="50%")
        ax[0, 1].scatter(xdate, q3, c="green", s=2, marker="o", label="75%")
        ax[0, 1].scatter(xdate, q4, c="black", s=2, marker="v", label="95%")
        ax[0, 1].set_ylabel(r"Box values, $\langle ADU \rangle_{norm}$")
        ax[0, 1].set_title("Amplifier boxes values vs nite")
        ax[0, 1].xaxis.set_major_locator(months)
        ax[0, 1].xaxis.set_major_formatter(dateFmt)
        ax[0, 1].xaxis.set_minor_locator(days)
        ax[0, 1].autoscale_view()
        kw_grid1 = {
            "color": "lightgray",
            "linestyle": "dotted",
            "dash_capstyle": "round",
            "alpha": 0.7,}
        ax[0, 1].grid(**kw_grid1)
        xlabels_ax1 = ax[0, 1].get_xticklabels()
        plt.setp(xlabels_ax1, rotation=30, fontsize=10)
        # Setup legend
        handles_ax1, labels_ax1 = ax[0, 1].get_legend_handles_labels()
        kw_lab1 = {
            "loc": "lower center",
            "ncol": len(labels_ax1),
            "markerscale": 2.1,
            "frameon": True,
            "fancybox": True,
            "framealpha": 0.5,
            "fontsize": 9,}
        ax[0, 1].legend(handles_ax1, labels_ax1, **kw_lab1)
        # For CCD3, AmpB, ylim [0.812, 0.820]
        # For CCD3, AmpA, ylim [1.0025, 1.0175]
        # For CCD2, AmpA, ylim [-0.00118, -0.0006]
        # For CCD3, AmpB, ylim [0.95, 1.2]
        # For CCD3, AmpA, ylim [0.96, 1.05]
        # For CCD3, AmpA, ylim [0.85, 1.3]
        # For CCD2, AmpA, ylim [-6, 65]
        # For CCD3, AmpB, ylim [0.95, 1.015]
        # ax[0, 1].set_ylim([-6, 65])
        #
        # 3rd density map (small scale when needed)
        rx0 = row_plot.compressed().ravel().min() 
        rx1 = row_plot.compressed().ravel().max()
        # For CCD3, ampB: ry0 = 0.812, ry1 = 0.820
        # For CCD3, ampA, ry0 = 1.0025, ry1 = 1.0175
        # For CCD3, ampA, ry0 = -0.00098, ry1 = -0.00088
        # For CCD3, AmpB, ry0 = 0.925, ry1 = 1.05
        # For CCD3, AmpA, ry0 = 0.98, ry1 = 1.02
        # For CCD3, AmpA, ry0 = 0.965, ry1 = 1.05
        # For CCD2, AmpA, ry0 = 1.5, ry1 = 14.3 
        # For CCD3, AmpB, ry0 = 0.95, ry1 = 1.04
        # For CCD2, AmpB, ry0 = -60000, ry1 = -2000
        ry0 = -45000
        ry1 = -2000
        kw_ax2 = {
            "bins": 80,
            "cmap": "viridis",
            "cmin": 1,
            "range": [[rx0, rx1], [ry0, ry1]],
            }
            # "range": [[0, 4096], [ry0, ry1]],
        ax[1, 0].axvline(x=row_issue1 * 16, c="orange")
        ax[1, 0].axvline(x=row_issue2 * 16, c="orange")
        p2 = ax[1, 0].hist2d(row_plot.compressed().ravel(), 
                            stt_sub["avg_n"].compressed().ravel(),
                            **kw_ax2)  
        ax[1, 0].set_xlabel("Row (box) number")
        ax[1, 0].set_ylabel(r"Boxes values, $\langle ADU \rangle_{norm}$")
        ax[1, 0].set_title("Density map for row profile".format(ccd))
        cb2 = plt.colorbar(p2[-1], ax=ax[1, 0])
        cb2.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb2.set_label("N")
        majorLocator_x2 = MultipleLocator(512)
        majorFormatter_x2 = FormatStrFormatter("%d")
        minorLocator_x2 = MultipleLocator(128)
        ax[1, 0].xaxis.set_minor_locator(minorLocator_x2)
        ax[1, 0].xaxis.set_major_locator(majorLocator_x2)
        ax[1, 0].xaxis.set_major_formatter(majorFormatter_x2)
        ax[1, 0].set_xlim([0,NBOX])
        #
        # 4th boxes for the issue region versus nite
        months = MonthLocator()
        days = DayLocator()
        dateFmt = DateFormatter("%b/%d/%y")
        months2 = MonthLocator()
        days2 = DayLocator()
        dateFmt2 = DateFormatter("%b/%d/%y")
        # Sort by expnum and transform nite to a plottable quantity
        xdate3 = [datetime.datetime.strptime(str(date[0]), "%Y%m%d")
                for date in np.sort(self.db["nite"])]
        kw_ax3 = {
            "s": 1,
            "c": "lightsteelblue",
            "marker": ".",
            }
        q0, q1, q2, q3, q4 = [], [], [], [], []
        # For each nite, plot profile
        for n in xrange(len(xdate3)):
            # x.compressed() is equivalent to x[~x.mask]
            tmp_stt = stt_sub["avg_n"][n,:].compressed()
            # To select only the issue region, must consider that one box
            # has been removed from the edge, and thus only 7 boxes remains
            tmp_stt = tmp_stt[row_issue1 * 7: row_issue2*7 + 1]
            tmp_date = [xdate3[n]] * tmp_stt.shape[0]
            q0.append(np.percentile(tmp_stt, 5)) 
            q1.append(np.percentile(tmp_stt, 25)) 
            q2.append(np.percentile(tmp_stt, 50))
            q3.append(np.percentile(tmp_stt, 75))
            q4.append(np.percentile(tmp_stt, 95)) 
            ax[1, 1].scatter(tmp_date, tmp_stt, **kw_ax3)
        ax[1, 1].scatter(xdate3, q0, c="violet", s=2, marker="^", label="5%")
        ax[1, 1].scatter(xdate3, q1, c="red", s=2, marker="s", label="25%")
        ax[1, 1].scatter(xdate3, q2, c="blue", s=2, marker="D", label="50%")
        ax[1, 1].scatter(xdate3, q3, c="green", s=2, marker="o", label="75%")
        ax[1, 1].scatter(xdate3, q4, c="black", s=2, marker="v", label="95%")
        ax[1, 1].set_ylabel(r"Box values, $\langle ADU \rangle_{norm}$")
        ax[1, 1].set_title("Issue region boxes values versus nite")
        ax[1, 1].xaxis.set_major_locator(months)
        ax[1, 1].xaxis.set_major_formatter(dateFmt)
        ax[1, 1].xaxis.set_minor_locator(days)
        ax[1, 1].autoscale_view()
        kw_grid3 = {
            "color": "lightgray",
            "linestyle": "dotted",
            "dash_capstyle": "round",
            "alpha": 0.7,}
        ax[1, 1].grid(**kw_grid3)
        # Setup xtick labels 
        xlabels_ax3 = ax[1, 1].get_xticklabels()
        plt.setp(xlabels_ax3, rotation=30, fontsize=10)
        # Setup legend
        handles_ax3, labels_ax3 = ax[1, 1].get_legend_handles_labels()
        kw_lab3 = {
            "loc": "lower center",
            "ncol": len(labels_ax3),
            "markerscale": 2.1,
            "frameon": True,
            "fancybox": True,
            "framealpha": 0.5,
            "fontsize": 9,}
        ax[1, 1].legend(handles_ax3, labels_ax3, **kw_lab3)
        # For CCD3, AmpB ylim [0.812, 0.820]
        # For CCD3, AmpA ylim [1.0025, 1.0175]
        # For CCD2, AmpA ylim [-0.001, -0.00085]
        # For CCD2, AmpB ylim [0.96, 1.2]
        # For CCD3, AmpA ylim [0.98, 1.03]
        # For CCD3, AmpA ylim [0.963, 1.1]
        # For CCD2, AmpA ylim [-2.5, 57.5]
        # For CCD3, AmpB ylim [0.95, 1.05]
        # ax[1, 1].set_ylim([-2.5, 57.5])
        # Final 
        supt = "[ZOOM] CCD{0}, amplifier showing".format(ccd)
        supt += " the 'specter'."
        supt += " Boxes of h=16 pix by w=128 pix"
        supt += "\nNOTE: Raw exposures with subtracted overscan."
        plt.suptitle(supt, color="blue")
        plt.subplots_adjust(left=0.1, bottom=0.08, top=0.9, right=0.98, 
                        hspace=0.26, wspace=0.26)
        if suffix is None:
            suffix = "pid" + str(os.getpid())
        if saveplot:
            outname = "rowProfileZOOM_{0:02}_{1}.png".format(ccd,suffix)
            plt.savefig(outname, dpi=400, facecolor="w", edgecolor="w", 
                    orientation="portrait", papertype=None, format="png", 
                    transparent=False, bbox_inches=None, pad_inches=0.1, 
                    frameon=None)
        plt.show()
    
    def histo_param(self, ccd, suffix=None, saveplot=True, param="norm"):
        """ Method to make histograms of the norm value, which was used to 
        normalize the boxes on which each CCD was divided. Or for other
        parameter
        Inputs:
        - ccd: ccd number to be plotted (0-99)
        - suffix: suffix to be append to the output filename
        - saveplot: wether to save or not an PNG output file
        - param: parameter to be plotted
        """
        print "Statistics array columns: {0}".format(self.stt.dtype.names)
        print "DB array columns: {0}".format(self.db.dtype.names)
        plt.close("all")
        fig, ax = plt.subplots(1, figsize=(8, 5))
        # Use the first element of each subset, because the rest are duplicates
        xnorm = np.copy(self.stt[param][:,0])
        aux_lab = "median of the CCD"
        kwh = {
            "bins": "doane",# "auto",
            "histtype": "stepfilled",
            "align": "mid",
            "orientation": "vertical",
            "log": False,
            "color": "royalblue",
            "label": aux_lab,
            }
        hist = ax.hist(xnorm, **kwh)
        colors = ["gray", "gold", "orange", "red", "gray"]
        for ind, qx in enumerate([5, 25, 50, 75, 95]):
            ax.axvline(x=np.percentile(xnorm, qx), c=colors[ind], 
                    label="{0}%".format(qx))
        if 1:
            # Ticks
            majorLocator_x = MultipleLocator(25)
            majorFormatter_x = FormatStrFormatter("%d")
            minorLocator_x = MultipleLocator(5)
            ax.xaxis.set_minor_locator(minorLocator_x)
            ax.xaxis.set_major_locator(majorLocator_x)
            ax.xaxis.set_major_formatter(majorFormatter_x)
        # Setup legend
        handles_ax, labels_ax = ax.get_legend_handles_labels()
        kw_lab = {
            "loc": "upper left",
            "ncol": len(labels_ax),
            "markerscale": 2.1,
            "frameon": True,
            "fancybox": True,
            "framealpha": 0.5,
            "fontsize": 9,}
        ax.legend(handles_ax, labels_ax, **kw_lab)
        # Labels and title
        ax.set_xlabel("Normalization value")
        ax.set_ylabel("N")
        title_aux = "Histogram for CCD{0:02} normalization value.".format(ccd)
        title_aux += " Note: raw w/o overscan subtraction"
        ax.set_title(title_aux)
        plt.subplots_adjust(left=0.07, bottom=0.1, right=0.98, top=0.93)
        if suffix is None:
            suffix = "pid" + str(os.getpid())
        if saveplot:
            outname = "normHisto_{0:02}_{1}.png".format(ccd,suffix)
            plt.savefig(outname, dpi=400, facecolor="w", edgecolor="w", 
                    orientation="portrait", papertype=None, format="png", 
                    transparent=False, bbox_inches=None, pad_inches=0.1, 
                    frameon=None)
        plt.show()
        

class One():
    def one_ccd(self, fnm):
        """Method to plot in a Log scale the CCD
        """
        from matplotlib.colors import LogNorm
        fits = FPBinned(fnm)
        for i in xrange(fits.ccd.shape[0]):
            plt.close("all")
            plt.imshow(fits.ccd[i, :, :], origin="lower", cmap="gray_r", 
                    norm=LogNorm())
            avg_region = np.mean(fits.ccd[i, 2208:3248, 100:1024])
            med_ccd = np.median(fits.ccd[i, :, :])
            perc = avg_region * 100. / med_ccd
            print "\n\t<region>={0}\n\tmedian(ccd)={1}\n\t{2}%".format(
                avg_region, med_ccd, perc)
            plt.show()


class Fit():
    def __init__(self):
        self.y = np.load("sub_adu.npy")
        self.x = np.load("sub_row.npy")


if __name__ == "__main__":
    print socket.gethostname()
    
    # Issue region  
    # CCD    |    Amplifier
    # ---------------------
    # 2      |    left (B)
    # 3      |    left (B)
    # 1      |    not victim
    # 6      |    not victim
    CCD = 3
    
    # Histo = Tide(stat="stat_{0:02}_16x128_y4e1_medN.npy".format(CCD),
    #         db="info_{0:02}_16x128_y4e1_medN.npy".format(CCD))
    # Histo = Tide(stat="stat_{0:02}_xtalked_16x128_medN.npy".format(CCD),
    #         db="info_{0:02}_xtalked_16x128_medN.npy".format(CCD))
    Histo = Tide(stat="stat_{0:02}_noOsc_16x128_medN.npy".format(CCD),
             db="info_{0:02}_noOsc_16x128_medN.npy".format(CCD))
    Histo.histo_param(ccd=CCD, suffix="noOsc_16x128_medN")

    # For RAW exposures with subtracted overscan and classical parameters
    # Tii = Tide(stat="stat_{0:02}_xtalked_16x128_medN.npy".format(CCD),
    #         db="info_{0:02}_xtalked_16x128_medN.npy".format(CCD))
    # Tii.zoom_row_profile(ccd=CCD, amp=0, suffix="xtalked_16x128_medN")

    # For REDUCED exposures, used for construct Y4E1 skytemplate
    # Ti = Tide(stat="stat_{0:02}_16x128_y4e1_medN.npy".format(CCD),
    #         db="info_{0:02}_16x128_y4e1_medN.npy".format(CCD)) 
    # 
    # Ti.zoom_row_profile(ccd=CCD, amp=1, suffix="ampA_16x128_y4e1_medN")
    # Ti.zoom_row_profile(ccd=CCD, amp=0, suffix="16x128_y4e1_medN")
    # Ti.row_profile(ccd=CCD, amp=0, suffix="16x128_y4e1_medN")
    #
    # For RAW exposures without overscan subtraction 
    # T = Tide(stat="stat_{0:02}_noOsc_16x128_medN.npy".format(CCD), 
    #         db="info_{0:02}_noOsc_16x128_medN.npy".format(CCD))        
    #
    # T.row_profile(ccd=CCD, suffix="NoOversc_16x128_medN")
    # T.zoom_row_profile(ccd=CCD, suffix="NoOversc_16x128_medN")
    #
    # T.row_profile(ccd=CCD, amp=1, suffix="NoOversc_AmpA_medN")
    # T.zoom_row_profile(ccd=CCD, amp=1, suffix="NoOversc_AmpA_medN")
    
    # One().one_ccd("Y4A1_20160801t1215_g_c03_r2930p01_skypca-tmpl.fits")

