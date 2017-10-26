"""Modifiable script to time series plotting, using 2D hitogram to trace
overdensities.
Customized fot plotting specific results, not fro general purposes.
"""
import sys
import os
import gc
import socket
import pandas as pd
import numpy as np
from numpy import ma
import datetime
import logging
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


class Compare():
    """ Ad-hoc class to plot side-by-side the profiles and results on the same
    scale, for the different reductions
    """

    def sel_amplifier(self, data, amp=None, ccdnum=None, NCCD=623, NBOX=4096,
                      column=None, Nrow=256, Ncol=16):
        """ Receives the entire CCD, and select only one amplifier. Returns
        a masked array of the original.
        Inputs
        - amp: str, either 'A' or 'B'
        - ccdnum: int. Both ccdnum and amp are used to get the AMPSEC 
        - NCCD: number of CCDs to mask. 
        - NBOX: number of boxes per CCD. Both NBOX and NCCD are used to 
        construct the mask
        - column: None, list or int. Selection of columns to be used for mask
        construction. Goues from 0 to Ncol-1
        - Nrow: Number of rows on which the CCD is divided
        - Ncol: Number of columns on which the CCD is divided
        """
        print "Using (expecting) {0} CCDs".format(NCCD)
        # Here use the table to get if the amplifier is at the left or right
        # in the CCD array
        if (amp.upper() in ["A", "B"]):
            amp = amp.upper()
            print "Amplifier {0} selected".format(amp)
        else:
            logging.error("Must input either \'A\' or \'B\' for the amplifier")
            exit(1)
        # Load the CCD sections table
        df_ccd = pd.read_csv("DECam_CCDsections.csv")
        # Using CCD number and AMPSECA, AMPSECB keywords
        if (amp == "A"):
            aux_amp = df_ccd.loc[(df_ccd["CCDNUM"] == ccdnum), "AMPSECA"]
            aux_amp = aux_amp.values[0]
        elif (amp == "B"):
            aux_amp = df_ccd.loc[(df_ccd["CCDNUM"] == ccdnum), "AMPSECB"]
            aux_amp = aux_amp.values[0]
        else:
            logging.error("Hmmm... something wrong happen")
        # Word of caution:
        # The AMPSECA/AMPSECB besides the sections, also give us the 
        # direction on which the data reading occurs. This is already taken 
        # into account whe the XTALK is performed. Look at this to better
        # understand the data
        # Now, the first dimension will give us the information about which
        # amplifier is to the left and which is to the right. Remember CCD
        # data dimensions is 2048 x 4096
        X, Y = aux_amp.strip("[").strip("]").split(",")
        X, Y = map(int, X.split(":")), map(int, Y.split(":"))
        if (max(X) == 1024):
            amp_lr = 0
        elif (max(X) == 2048):
            amp_lr = 1
        else:
            logging.error("Error when deciding AMP order in CCD")
            exit(1)
        # Distinguish per amplifier
        # Giving the binnig of 128 x 16 pixels, the amount of rows (larger dim)
        # is 256, and the amount of columns (shorter dim) is 16
        row_ini, row_end = 0, Nrow
        # By eye, we can say the bump is mostly located at rows = 150 to 202
        if amp_lr == 0:
            # left amplifier without the edge box
            left = np.arange(row_ini * Ncol, row_end * Ncol)[::Ncol]
            idx_left = []
            for L in left:
                if (column is None):
                    # avoid the edge
                    idx_left += list(np.arange(L + 1, L + int(Ncol/2)))
                elif isinstance(column, int):
                    # select only one column
                    idx_left += [L + column]
                elif isinstance(column, list):
                    for item in column:
                        idx_left += [L + item]
                else:
                    logging.error("Column must be either None, list or int")
                    exit(1)
            idx = np.array(idx_left)
            # sort the array
            idx = np.sort(idx)
        elif amp_lr == 1:
            # right amplifier without the edge box
            right = np.arange(row_ini*Ncol + int(Ncol/2), 
                              row_end*Ncol + int(Ncol/2))[::Ncol]
            idx_right = []
            for R in right:
                if (column is None):
                    # avoid the edge
                    idx_right += list(np.arange(R, R + int(Ncol/2)))
                elif isinstance(column, int):
                    # select only one column
                    idx_right += [R + column]
                elif isinstance(column, list):
                    for item in column:
                        idx_right += [R + item]
                else:
                    logging.error("Column must be either None, list or int")
                    exit(1)
            idx = np.array(idx_right)
            idx = np.sort(idx)
        else:
            logging.error("Error in amplifier selection")
            exit(1)
        
        # Define the masks using the above indices
        row_mask_aux = np.ones((NCCD, NBOX), dtype=bool)
        row_mask_aux[:, idx] = False
        row_plot = np.tile(np.arange(NBOX), NCCD).reshape(NCCD, NBOX)
        row_plot = ma.array(row_plot, mask=row_mask_aux)
        # For the boxes, I need to select only left side, using a mask
        # Apply the mask to the boxes-statistics and db-info arrays
        not_idx = np.ones((NCCD, NBOX),dtype=bool)
        not_idx[:, idx] = False
        data_sub = ma.array(data, mask=not_idx)
        return data_sub, row_plot

    def side1(self, ccd=None, amp=0, saveplot=True, suffix=None):
        """ Options of plot 1, side by side per CCD. All at the same scale
        or at least compatible by some shown factor (put some annotation)
        Upper row contains the profile for the ccd rows (larger dimension),
        zoomed to the region of higher density of points
        Lower row contains the histogram of the values used to normalize the
        boxes, usually the median of the values per CCD.
        """
        # For NoOsc
        # fnm1 = "stat_{0:02}_noOsc_16x128_medN.npy".format(CCD)
        # For ccd02NoOsc
        fnm1 = "stat_{0:02}_ccd02NoOsc_medN.npy".format(CCD)
        fnm2 = "stat_{0:02}_xtalked_16x128_medN.npy".format(CCD)
        fnm3 = "stat_{0:02}_16x128_y4e1_medN.npy".format(CCD)
        fnm = [fnm1, fnm2, fnm3]
        # To create variable names from string
        this_module = sys.modules[__name__]
        for index,item in enumerate(fnm):
            name = "aux_test_{0}".format(index + 1)
            setattr(this_module, name, index)
        # Retrieve only one amplifier, masked arrays
        tmp1, tmp2, tmp3 = np.load(fnm1), np.load(fnm2), np.load(fnm3)
        row_plot = None
        d1, row_plot1 = Compare().sel_amplifier(tmp1, amp=amp, ccdnum=ccd)
        d2, row_plot2 = Compare().sel_amplifier(tmp2, amp=amp, ccdnum=ccd)
        d3, row_plot3 = Compare().sel_amplifier(tmp3, amp=amp, ccdnum=ccd)
        # Define plot
        plt.close("all")
        fig, ax = plt.subplots(2, 3, figsize=(12, 8))
        # Density histograms
        # For CCD3, ywide 0.1
        ywide = .1
        kw0 = {
            "bins": 100,
            "cmap": "viridis",
            "cmin": 1,
            }
        # For CCD3 ampB
        # kw0.update({"range": [[0, 4096] , [0.807, 0.807 + ywide/4.]],})
        # For CCD3 ampA
        # kw0.update({"range": [[0, 4096] , [0.995, 0.995 + ywide/4.]],})
        # For CCD3 ampB, ccd02NoOsc
        # kw0.update({"range": [[0, 4096] , [0.925, 0.925 + ywide]],})
        kw0.update({"range": [[0, 4096] , [0.96, 0.96 + ywide]],})
        p1 = ax[0, 0].hist2d(row_plot.compressed().ravel(),
                            d1["avg_n"].compressed().ravel(),
                            **kw0)
        # For CCD3 ampB
        # kw0.update({"range": [[0, 4096] , [0.94, 0.94 + ywide]],})
        # For ccd3 ampA
        # kw0.update({"range": [[0, 4096] , [0.95, 0.95 + ywide]],})
        kw0.update({"range": [[0, 4096] , [0.95, 0.95 + ywide]],})
        p2 = ax[0, 1].hist2d(row_plot.compressed().ravel(),
                            d2["avg_n"].compressed().ravel(),
                            **kw0)
        # For CCD3 ampB
        # kw0.update({"range": [[0, 4096] , [0.96, 0.96 + ywide]],})
        # For ccd3 ampA
        # kw0.update({"range": [[0, 4096] , [0.95, 0.95 + ywide]],})
        kw0.update({"range": [[0, 4096] , [0.95, 0.95 + ywide]],})
        p3 = ax[0, 2].hist2d(row_plot.compressed().ravel(),
                            d3["avg_n"].compressed().ravel(),
                            **kw0)
        # Norm histograms
        parameter = "norm"
        xnorm1 = tmp1[parameter][:, 0]
        xnorm2 = tmp2[parameter][:, 0]
        xnorm3 = tmp3[parameter][:, 0]
        kw_hist = {
            "bins": "doane",# "auto",
            "histtype": "stepfilled",
            "align": "mid",
            "orientation": "vertical",
            "log": False,
            "color": "gray",
            }
        h1 = ax[1, 0].hist(xnorm1, **kw_hist)
        h2 = ax[1, 1].hist(xnorm2, **kw_hist)
        h3 = ax[1, 2].hist(xnorm3, **kw_hist)
        colors = ["gold", "lime", "cyan", "red", "gold"]
        for ind, qx in enumerate([5, 25, 50, 75, 95]):
            ax[1, 0].axvline(x=np.percentile(xnorm1, qx), c=colors[ind],
                            label="{0}%".format(qx))
            ax[1, 1].axvline(x=np.percentile(xnorm2, qx), c=colors[ind],
                            label="{0}%".format(qx))
            ax[1, 2].axvline(x=np.percentile(xnorm3, qx), c=colors[ind],
                            label="{0}%".format(qx))
        # Setup legend for histograms
        handles_ax1, labels_ax1 = ax[1, 0].get_legend_handles_labels()
        handles_ax2, labels_ax2 = ax[1, 1].get_legend_handles_labels()
        handles_ax3, labels_ax3 = ax[1, 2].get_legend_handles_labels()
        kw_lab = {
            "loc": "upper right",
            "ncol": 1,
            "markerscale": 2.1,
            "frameon": True,
            "fancybox": True,
            "framealpha": 0.5,
            "fontsize": 9,}
        ax[1, 0].legend(handles_ax1, labels_ax1, **kw_lab)
        ax[1, 1].legend(handles_ax2, labels_ax1, **kw_lab)
        ax[1, 2].legend(handles_ax3, labels_ax1, **kw_lab)
        # Colorbars
        kw_cbar = {"orientation": "horizontal", "aspect": 50}
        cb1 = plt.colorbar(p1[-1], ax=ax[0, 0], **kw_cbar)
        cb1.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb1.set_label("N")
        cb2 = plt.colorbar(p2[-1], ax=ax[0, 1], **kw_cbar)
        cb2.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb2.set_label("N")
        cb3 = plt.colorbar(p3[-1], ax=ax[0, 2], **kw_cbar)
        cb3.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb3.set_label("N")
        # Setup ticks
        majorLocator_x0 = MultipleLocator(1024)
        majorFormatter_x0 = FormatStrFormatter("%d")
        minorLocator_x0 = MultipleLocator(128)
        ax[0, 0].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 0].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 0].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 1].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 1].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 1].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 2].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 2].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 2].xaxis.set_major_formatter(majorFormatter_x0)
        # Labels
        kw_lab = {"fontsize": 10, "color": "blue",}
        ax[0, 0].set_title("xtalked, CCD2 null overscan. Scale: 1X", **kw_lab)
        ax[0, 1].set_title("xtalked. Scale: 1X", **kw_lab)
        ax[0, 2].set_title("reduced. Scale: 1X", **kw_lab)
        ax[0, 0].set_ylabel(r"$<x_{box}>_{norm}$")
        ax[0 ,0].set_xlabel("rows")
        ax[0, 1].set_xlabel("rows")
        ax[0, 2].set_xlabel("rows")
        ax[1, 0].set_ylabel("N")
        ax[1, 0].set_xlabel("Norm value")
        ax[1, 1].set_xlabel("Norm value")
        ax[1, 2].set_xlabel("Norm value")
        # For CCD3 plt.suptitle("CCD {0}, issue Amp".format(CCD))
        plt.suptitle("CCD {0}, *non* issue Amp".format(CCD))

        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.98, top=0.90,
                        wspace=0.16, hspace=0.07)
        if suffix is None:
            suffix = "pid" + str(os.getpid())
        if saveplot:
            outname = "profileComp_{0:02}_{1}.png".format(ccd,suffix)
            plt.savefig(outname, dpi=400, facecolor="w", edgecolor="w",
                    orientation="portrait", papertype=None, format="png",
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    frameon=None)
        plt.show()

    def side2(self, ccd=None, amp=0, saveplot=True, suffix=None):
        """ Options of plot 2, side by side per CCD, only 2 sets
        """
        # For CCD2 NoOsc
        # fnm1 = "stat_{0:02}_noOsc_16x128_medN.npy".format(CCD)
        fnm1 = "stat_{0:02}_ccd02NoOsc_medN.npy".format(CCD)
        fnm2 = "stat_{0:02}_xtalked_16x128_medN.npy".format(CCD)
        fnm = [fnm1, fnm2]
        # Retrieve only one amplifier, masked arrays
        tmp1, tmp2 = np.load(fnm1), np.load(fnm2)
        row_plot = None
        d1, row_plot1 = Compare().sel_amplifier(tmp1, amp=amp, ccdnum=ccd)
        d2, row_plot2 = Compare().sel_amplifier(tmp2, amp=amp, ccdnum=ccd)
        # Define plot
        plt.close("all")
        fig, ax = plt.subplots(2, 2, figsize=(8, 8))
        # Density histograms
        # For CCD3, ywide 0.1
        # For CCD2, ywide = 6400
        # For CCD2, ampA ccd02NoOsc ywide = 4000
        ywide = 7E3
        kw0 = {
            "bins": 100,
            "cmap": "viridis",
            "cmin": 1,
            }
        # For CCD3
        # kw0.update({"range": [[0, 4096] , [0.807, 0.807 + ywide/4.]],})
        # kw0.update({"range": [[0, 4096] , [0.995, 0.995 + ywide/4.]],})
        # For CCD2
        # kw0.update({"range": [[0, 4096] , [2.34, 2.34 + ywide/4000.]],})
        # kw0.update({"range": [[0, 4096] , [-0.001, -0.00064]],})
        # For CCD2, ampA, ccd02NoOsc
        # kw0.update({"range": [[0, 4096] , [-0.0011, -0.0011 + ywide]],})
        p1 = ax[0, 0].hist2d(row_plot.compressed().ravel(),
                            d1["avg_n"].compressed().ravel(),
                            **kw0)
        # For CCD3
        # kw0.update({"range": [[0, 4096] , [0.94, 0.94 + ywide]],})
        # For CCD2
        # kw0.update({"range": [[0, 4096] , [-8000, -8000 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [1.8, 5]],})
        # For CCD2, ampA, ccd02NoOsc
        # kw0.update({"range": [[0, 4096] , [1.8, 1.8 + 1E4 * ywide]],})
        kw0.update({"range": [[0, 4096] , [-9500, -9500 + ywide]],})
        p2 = ax[0, 1].hist2d(row_plot.compressed().ravel(),
                            d2["avg_n"].compressed().ravel(),
                            **kw0)
        # Norm histograms
        parameter = "norm"
        xnorm1 = tmp1[parameter][:, 0]
        xnorm2 = tmp2[parameter][:, 0]
        kw_hist = {
            "bins": "doane",# "auto",
            "histtype": "stepfilled",
            "align": "mid",
            "orientation": "vertical",
            "log": False,
            "color": "gray",
            }
        h1 = ax[1, 0].hist(xnorm1, **kw_hist)
        h2 = ax[1, 1].hist(xnorm2, **kw_hist)
        colors = ["gold", "lime", "cyan", "red", "gold"]
        for ind, qx in enumerate([5, 25, 50, 75, 95]):
            ax[1, 0].axvline(x=np.percentile(xnorm1, qx), c=colors[ind],
                            label="{0}%".format(qx))
            ax[1, 1].axvline(x=np.percentile(xnorm2, qx), c=colors[ind],
                            label="{0}%".format(qx))
        # Setup legend for histograms
        handles_ax1, labels_ax1 = ax[1, 0].get_legend_handles_labels()
        handles_ax2, labels_ax2 = ax[1, 1].get_legend_handles_labels()
        kw_lab = {
            "loc": "upper left",
            "ncol": 1,
            "markerscale": 2.1,
            "frameon": True,
            "fancybox": True,
            "framealpha": 0.5,
            "fontsize": 9,}
        ax[1, 0].legend(handles_ax1, labels_ax1, **kw_lab)
        ax[1, 1].legend(handles_ax2, labels_ax1, **kw_lab)
        # Colorbars
        kw_cbar = {"orientation": "horizontal", "aspect": 50}
        cb1 = plt.colorbar(p1[-1], ax=ax[0, 0], **kw_cbar)
        cb1.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb1.set_label("N")
        cb2 = plt.colorbar(p2[-1], ax=ax[0, 1], **kw_cbar)
        cb2.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb2.set_label("N")
        # Setup ticks
        majorLocator_x0 = MultipleLocator(1024)
        majorFormatter_x0 = FormatStrFormatter("%d")
        minorLocator_x0 = MultipleLocator(128)
        ax[0, 0].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 0].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 0].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 1].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 1].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 1].xaxis.set_major_formatter(majorFormatter_x0)
        # Labels
        kw_lab = {"fontsize": 10, "color": "blue",}
        # ax[0, 0].set_title("xtalked, no overscan", **kw_lab)
        ax[0, 0].set_title("xtalked, CCD2 null overscan. Scale: 5,200X",
                           **kw_lab)
        ax[0, 1].set_title("xtalked. Scale: 1X", **kw_lab)
        ax[0, 0].set_ylabel(r"$<x_{box}>_{norm}$")
        ax[0 ,0].set_xlabel("rows")
        ax[0, 1].set_xlabel("rows")
        ax[1, 0].set_ylabel("N")
        ax[1, 0].set_xlabel("Norm value")
        ax[1, 1].set_xlabel("Norm value")
        # For CCD3 plt.suptitle("CCD {0}, issue Amp".format(CCD))
        plt.suptitle("CCD {0}, issue Amp".format(CCD))

        plt.subplots_adjust(left=.1, bottom=0.06, right=0.98, top=0.90,
                        wspace=0.28, hspace=0.07)
        if suffix is None:
            suffix = "pid" + str(os.getpid())
        if saveplot:
            outname = "profileComp_{0:02}_{1}.png".format(ccd,suffix)
            plt.savefig(outname, dpi=400, facecolor="w", edgecolor="w",
                    orientation="portrait", papertype=None, format="png",
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    frameon=None)
        plt.show()

    def side3(self, ccd=None, amp=0, saveplot=True, suffix=None):
        """ Options of plot settings #3, 4 columns, side by side per CCD.
        All at the same scale
        or at least compatible by some shown factor (put some annotation)
        Upper row contains the profile for the ccd rows (larger dimension),
        zoomed to the region of higher density of points
        Lower row contains the histogram of the values used to normalize the
        boxes, usually the median of the values per CCD.
        """
        # For NoOsc
        # fnm1 = "stat_{0:02}_noOsc_16x128_medN.npy".format(CCD)
        # For ccd02NoOsc
        # fnm1 = "stat_{0:02}_noOsc02_16x128_notN_v1.npy".format(CCD)
        # fnm2 = "stat_{0:02}_noOsc02_16x128_notN_v2.npy".format(CCD)
        # fnm3 = "stat_{0:02}_noOsc02_16x128_notN_v3.npy".format(CCD)
        fnm1 = "stat_{0:02}_16x128_medN_v2.npy".format(CCD)
        fnm2 = "stat_{0:02}_16x128_medN_v2.A.npy".format(CCD)
        fnm3 = "stat_{0:02}_16x128_medN_v2.B.npy".format(CCD)
        fnm4 = "stat_{0:02}_16x128_medN_normal.npy".format(CCD)
        # fnm1 = "stat_{0:02}_16x128_notN_v2.npy".format(CCD)
        # fnm2 = "stat_{0:02}_16x128_notN_v2.A.npy".format(CCD)
        # fnm3 = "stat_{0:02}_16x128_notN_v2.B.npy".format(CCD)
        # fnm4 = "stat_{0:02}_16x128_notN_normal.npy".format(CCD)
        fnm = [fnm1, fnm2, fnm3, fnm4]
        # To create variable names from string
        this_module = sys.modules[__name__]
        for index,item in enumerate(fnm):
            name = "aux_test_{0}".format(index + 1)
            setattr(this_module, name, index)
        print "{0}".format(this_module)
        # Retrieve only one amplifier, masked arrays
        tmp1, tmp2, tmp3 = np.load(fnm1), np.load(fnm2), np.load(fnm3)
        tmp4 = np.load(fnm4)
        row_plot = None
        d1, row_plot1 = Compare().sel_amplifier(tmp1, amp=amp, ccdnum=ccd)
        d2, row_plot2 = Compare().sel_amplifier(tmp2, amp=amp, ccdnum=ccd)
        d3, row_plot3 = Compare().sel_amplifier(tmp3, amp=amp, ccdnum=ccd)
        d4, row_plot4 = Compare().sel_amplifier(tmp4, amp=amp, ccdnum=ccd)
        # Define plot
        plt.close("all")
        fig, ax = plt.subplots(2, 5, figsize=(15, 8))
        # Density histograms
        # For CCD3, ywide 0.1
        ywide = .1
        # ywide = 60
        kw0 = {
            "bins": 100,
            "cmap": "viridis",
            "cmin": 1,
            }
        yvalue = "med_n"
        # NOTE: to plot the med_n is OK even when no normalization is selected
        #       from the initial .npy file, because when norm=none then norm=1
        #
        # For CCD3 ampB
        # kw0.update({"range": [[0, 4096] , [0.807, 0.807 + ywide/4.]],})
        # For CCD3 ampA
        # kw0.update({"range": [[0, 4096] , [0.995, 0.995 + ywide/4.]],})
        # For CCD3 ampB, ccd02NoOsc
        # kw0.update({"range": [[0, 4096] , [0.925, 0.925 + ywide]],})
        kw0.update({"range": [[0, 4096] , [0.95, 0.95 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [60, 60 + ywide]],})
        p1 = ax[0, 0].hist2d(row_plot1.compressed().ravel(),
                             d1[yvalue].compressed().ravel(),
                             **kw0)
        # For CCD3 ampB
        # kw0.update({"range": [[0, 4096] , [0.94, 0.94 + ywide]],})
        # For ccd3 ampA
        # kw0.update({"range": [[0, 4096] , [0.95, 0.95 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [0.95, 0.95 + ywide]],})
        p2 = ax[0, 1].hist2d(row_plot2.compressed().ravel(),
                             d2[yvalue].compressed().ravel(),
                             **kw0)
        # For CCD3 ampB
        # kw0.update({"range": [[0, 4096] , [0.96, 0.96 + ywide]],})
        # For ccd3 ampA
        # kw0.update({"range": [[0, 4096] , [0.95, 0.95 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [0.95, 0.95 + ywide]],})
        p3 = ax[0, 2].hist2d(row_plot3.compressed().ravel(),
                             d3[yvalue].compressed().ravel(),
                             **kw0)
        # Normal processing
        p4 = ax[0, 3].hist2d(row_plot4.compressed().ravel(),
                             d4[yvalue].compressed().ravel(),
                             **kw0)
        # Difference between normal and some coefficient
        # kw0.update({"range": [[0, 4096] , [-200, 200]],})
        kw0.update({"range": [[0, 4096] , [-5, 5]],})
        deltax = row_plot3.compressed().ravel()
        deltay = (d4[yvalue].compressed() - d3[yvalue].compressed()).ravel()
        # deltay = (d4[yvalue] - d2[yvalue]).compressed().ravel()
        kw0.update({"cmap": "plasma"})
        p5 = ax[0, 4].hist2d(deltax, deltay, **kw0)
        # Some statistics histogram
        # dtype=[('norm', '<f4'), ('med', '<f4'), ('avg', '<f4'),
        # ('med_n', '<f4'), ('avg_n', '<f4'), ('rms_n', '<f4'),
        # ('unc_n', '<f4'), ('mad_n', '<f4')])
        parameter = "med_n"
        xnorm1 = tmp1[parameter][:, 0]
        xnorm2 = tmp2[parameter][:, 0]
        xnorm3 = tmp3[parameter][:, 0]
        xnorm4 = tmp4[parameter][:, 0]
        kw_hist = {
            "bins": "auto",#""doane",# "auto",
            "histtype": "stepfilled",
            "align": "mid",
            "orientation": "vertical",
            "log": False,
            "color": "gray",
            }
        h1 = ax[1, 0].hist(xnorm1, **kw_hist)
        h2 = ax[1, 1].hist(xnorm2, **kw_hist)
        h3 = ax[1, 2].hist(xnorm3, **kw_hist)
        h4 = ax[1, 3].hist(xnorm4, **kw_hist)
        colors = ["gold", "lime", "cyan", "red", "gold"]
        for ind, qx in enumerate([5, 25, 50, 75, 95]):
            ax[1, 0].axvline(x=np.percentile(xnorm1, qx), c=colors[ind],
                            label="{0}%".format(qx))
            ax[1, 1].axvline(x=np.percentile(xnorm2, qx), c=colors[ind],
                            label="{0}%".format(qx))
            ax[1, 2].axvline(x=np.percentile(xnorm3, qx), c=colors[ind],
                            label="{0}%".format(qx))
            ax[1, 3].axvline(x=np.percentile(xnorm4, qx), c=colors[ind],
                            label="{0}%".format(qx))
        # For the separate histogram
        # If plotting the norm, it is common for the entire CCD
        # If plotting other parameter, keep in mind to change setup
        # parameter = "norm"
        # xnorm5 = tmp4[parameter][:, 0]
        kw_hist = {
            "bins": "auto",#""doane",# "auto",
            "histtype": "stepfilled",
            "align": "mid",
            "orientation": "vertical",
            "log": False,
            "color": "gray",
            }
        h5 = ax[1, 4].hist(deltay, **kw_hist)
        colors = ["gold", "lime", "cyan", "red", "gold"]
        for ind, qx in enumerate([5, 25, 50, 75, 95]):
            ax[1, 4].axvline(x=np.percentile(deltay, qx), c=colors[ind],
                            label="{0}%".format(qx))
        # Setup plot limits for histograms
        # ax[1, 0].set_xlim([60, 200])
        # ax[1, 1].set_xlim([60, 200])
        # ax[1, 2].set_xlim([60, 200])
        # ax[1, 3].set_xlim([60, 200])
        # ax[1, 4].set_xlim([-100, 100])
        ax[1, 0].set_xlim([1.2, 1.5])
        ax[1, 1].set_xlim([1.2, 1.5])
        ax[1, 2].set_xlim([1.2, 1.5])
        ax[1, 3].set_xlim([1.2, 1.5])
        ax[1, 4].set_xlim([-0.1, 0.1])
        # ax[1, 0].set_xlim([70, 220])
        # ax[1, 1].set_xlim([70, 220])
        # ax[1, 2].set_xlim([70, 220])
        # Setup legend for histograms
        handles_ax1, labels_ax1 = ax[1, 0].get_legend_handles_labels()
        handles_ax2, labels_ax2 = ax[1, 1].get_legend_handles_labels()
        handles_ax3, labels_ax3 = ax[1, 2].get_legend_handles_labels()
        handles_ax4, labels_ax4 = ax[1, 3].get_legend_handles_labels()
        handles_ax5, labels_ax5 = ax[1, 3].get_legend_handles_labels()
        kw_lab = {
            "loc": "upper right",
            "ncol": 1,
            "markerscale": 2.1,
            "frameon": True,
            "fancybox": True,
            "framealpha": 0.5,
            "fontsize": 9,}
        ax[1, 0].legend(handles_ax1, labels_ax1, **kw_lab)
        ax[1, 1].legend(handles_ax2, labels_ax2, **kw_lab)
        ax[1, 2].legend(handles_ax3, labels_ax3, **kw_lab)
        ax[1, 3].legend(handles_ax4, labels_ax4, **kw_lab)
        ax[1, 4].legend(handles_ax5, labels_ax5, **kw_lab)
        # Colorbars
        kw_cbar = {"orientation": "horizontal", "aspect": 50}
        cb1 = plt.colorbar(p1[-1], ax=ax[0, 0], **kw_cbar)
        cb1.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb1.set_label("N")
        cb2 = plt.colorbar(p2[-1], ax=ax[0, 1], **kw_cbar)
        cb2.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb2.set_label("N")
        cb3 = plt.colorbar(p3[-1], ax=ax[0, 2], **kw_cbar)
        cb3.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb3.set_label("N")
        cb4 = plt.colorbar(p4[-1], ax=ax[0, 3], **kw_cbar)
        cb4.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb4.set_label("N")
        cb5 = plt.colorbar(p5[-1], ax=ax[0, 4], **kw_cbar)
        cb5.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
        cb5.set_label("N")
        # Fontsize of the numbers!
        cb1.ax.tick_params(labelsize=7)
        cb2.ax.tick_params(labelsize=7)
        cb3.ax.tick_params(labelsize=7)
        cb4.ax.tick_params(labelsize=7)
        cb5.ax.tick_params(labelsize=7)
        # Setup ticks
        majorLocator_x0 = MultipleLocator(1024)
        majorFormatter_x0 = FormatStrFormatter("%d")
        minorLocator_x0 = MultipleLocator(128)
        ax[0, 0].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 0].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 0].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 1].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 1].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 1].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 2].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 2].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 2].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 3].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 3].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 3].xaxis.set_major_formatter(majorFormatter_x0)
        # For CCD3 plt.suptitle("CCD {0}, issue Amp".format(CCD))
        ttext = "CCD {0}, amp B.".format(CCD)
        ttext += " Comparison of xtalk coeffs variations. Overscan for all."
        ttext += "\nCCD-median normalized imgs. (victim-source)"
        plt.suptitle(ttext)
        # Labels
        kw_lab = {"fontsize": 10, "color": "blue",}
        ax[0, 0].set_title("vers2. (3B-2B)/20", **kw_lab)
        ax[0, 1].set_title("vers2.A zero coeff 3B-2B", **kw_lab)
        ax[0, 2].set_title("vers2.B (3B-2B)*20", **kw_lab)
        ax[0, 3].set_title("normal process", **kw_lab)
        kw_lab.update({"color": "red"})
        ax[0, 4].set_title("DIFF: (normal-20x)", **kw_lab)
        ax[0, 0].set_ylabel(r"$med({box})$")
        ax[0 ,0].set_xlabel("rows")
        ax[0, 1].set_xlabel("rows")
        ax[0, 2].set_xlabel("rows")
        ax[0, 3].set_xlabel("rows")
        ax[0, 4].set_xlabel("rows")
        ax[1, 0].set_ylabel("N")
        ax[1, 0].set_xlabel(r"$med({box})$")
        ax[1, 1].set_xlabel(r"$med({box})$")
        ax[1, 2].set_xlabel(r"$med({box})$")
        ax[1, 3].set_xlabel(r"$med({box})$")
        ax[1, 4].set_xlabel(r"DIFF: (normal-20x)")
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.98, top=0.90,
                        wspace=0.16, hspace=0.07)
        if suffix is None:
            suffix = "pid" + str(os.getpid())
        if saveplot:
            outname = "profileComp_{0:02}_{1}.png".format(ccd,suffix)
            plt.savefig(outname, dpi=400, facecolor="w", edgecolor="w",
                    orientation="portrait", papertype=None, format="png",
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    frameon=None)
        plt.show()
    
    def side4(self, ccd=None, amp=0, saveplot=True, suffix=None, region=None,
              column=None):
        """ Options of plot settings #4, 6 columns, side by side per CCD.
        All at the same scale
        or at least compatible by some shown factor (put some annotation)
        Upper row contains the profile for the ccd rows (larger dimension),
        zoomed to the region of higher density of points
        Lower row contains the histogram of the values used to normalize the
        boxes, usually the median of the values per CCD.
        """
        if False:
            if (region is None):
                fnm1 = "stat_{0:02}_16x128_medN_v2.npy".format(CCD)
                fnm2 = "stat_{0:02}_16x128_medN_v2.A.npy".format(CCD)
                fnm3 = "stat_{0:02}_16x128_medN_v2.B.npy".format(CCD)
                fnm4 = "stat_{0:02}_16x128_medN_v2.C.npy".format(CCD)
                fnm5 = "stat_{0:02}_16x128_medN_v2.D.npy".format(CCD)
                fnm6 = "stat_{0:02}_16x128_medN_normal.npy".format(CCD)
            else:
                fnm1 = "stat_{0:02}_16x128_medN_v2_{1}.npy".format(CCD, region)
                fnm2 = "stat_{0:02}_16x128_medN_v2.A_{1}.npy".format(CCD, region)
                fnm3 = "stat_{0:02}_16x128_medN_v2.B_{1}.npy".format(CCD, region)
                fnm4 = "stat_{0:02}_16x128_medN_v2.C_{1}.npy".format(CCD, region)
                fnm5 = "stat_{0:02}_16x128_medN_v2.D_{1}.npy".format(CCD, region)
                fnm6 = "stat_{0:02}_16x128_medN_normal.npy".format(CCD)
            fnm = [fnm1, fnm2, fnm3, fnm4, fnm5, fnm6]
        if False:
            fnm1 = "stat_{0:02}_noOsc02_16x128_medN_v2.E_reg3.npy".format(CCD)
            fnm2 = "stat_{0:02}_noOsc02_16x128_medN_v2.F_reg3.npy".format(CCD)
            fnm3 = "stat_{0:02}_noOsc02_16x128_medN_v2.G_reg3.npy".format(CCD)
            fnm4 = "stat_{0:02}_noOsc02_16x128_medN_v2.H_reg3.npy".format(CCD)
            fnm5 = "stat_{0:02}_noOsc02_16x128_medN_v2.B_reg3.npy".format(CCD)
            # v1 has unaltered coefficients for ccd03B - ccd02B
            # fnm6 = "stat_{0:02}_noOsc02_16x128_medN_v1.npy".format(CCD)
            fnm6 = "stat_{0:02}_noOsc02_16x128_medN_unalt_reg3.npy".format(CCD)
        if False:
            fnm1 = "stat_{0:02}_noOsc02_16x128_medN_v2.E.npy".format(CCD)
            fnm2 = "stat_{0:02}_noOsc02_16x128_medN_v2.F.npy".format(CCD)
            fnm3 = "stat_{0:02}_noOsc02_16x128_medN_v2.G.npy".format(CCD)
            fnm4 = "stat_{0:02}_noOsc02_16x128_medN_v2.H.npy".format(CCD)
            fnm5 = "stat_{0:02}_noOsc02_16x128_medN_v2.B.npy".format(CCD)
            # v1 has unaltered coefficients for ccd03B - ccd02B
            fnm6 = "stat_{0:02}_noOsc02_16x128_medN_v1.npy".format(CCD)
        if False:
            fnm1 = "stat_{0:02}_noOsc02_16x128_medN_v2.G_{1}.npy".format(CCD, region)
            fnm2 = "stat_{0:02}_noOsc02_16x128_medN_v2.H_{1}.npy".format(CCD, region)
            fnm3 = "stat_{0:02}_noOsc02_16x128_medN_v2.B_{1}.npy".format(CCD, region)
            fnm4 = "stat_{0:02}_noOsc02_16x128_medN_v2.C_{1}.npy".format(CCD, region)
            fnm5 = "stat_{0:02}_noOsc02_16x128_medN_v2.D_{1}.npy".format(CCD, region)
            fnm6 = "stat_{0:02}_noOsc02_16x128_medN_unalt_{1}.npy".format(CCD, region) 
        fnm = [fnm1, fnm2, fnm3, fnm4, fnm5, fnm6]
        # To create variable names from string
        this_module = sys.modules[__name__]
        for index,item in enumerate(fnm):
            name = "aux_test_{0}".format(index + 1)
            setattr(this_module, name, index)
        print "{0}".format(this_module)
        # Retrieve only one amplifier, masked arrays
        tmp1, tmp2, tmp3 = np.load(fnm1), np.load(fnm2), np.load(fnm3)
        tmp4, tmp5 = np.load(fnm4), np.load(fnm5)
        tmp6 = np.load(fnm6)
        row_plot = None
        #
        # Here use the selection of column
        dsel = {"amp": amp, "ccdnum": ccd, "column":column}
        d1, row_plot1 = Compare().sel_amplifier(tmp1, **dsel)
        d2, row_plot2 = Compare().sel_amplifier(tmp2, **dsel)
        d3, row_plot3 = Compare().sel_amplifier(tmp3, **dsel)
        d4, row_plot4 = Compare().sel_amplifier(tmp4, **dsel)
        d5, row_plot5 = Compare().sel_amplifier(tmp5, **dsel)
        d6, row_plot6 = Compare().sel_amplifier(tmp6, **dsel)
        #
        # Define plot
        plt.close("all")
        fig, ax = plt.subplots(2, 6, figsize=(17, 8))
        # Density histograms
        # For CCD3, ywide 0.1
        # ywide = .1
        ywide = 0.1
        # ywide = 60
        kw0 = {
            "bins": 100,
            "cmap": "viridis", # "plasma", # "viridis",
            "cmin": 1,
            }
        yvalue = "med_n"
        # NOTE: to plot the med_n is OK even when no normalization is selected
        #       from the initial .npy file, because when norm=none then norm=1
        #
        # kw0.update({"range": [[0, 4096] , [0.95, 0.95 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [0.92, 0.92 + ywide]],})
        kw0.update({"range": [[0, 4096] , [0.98, 0.98 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [0.7, 0.7 + 3.*ywide]],})
        # Use an UNIQUE range of values for the colorbar
        kw0.update({"vmin": 0, "vmax": 250})
        # kw0.update({"vmin": 0, "vmax": 900})
        p1 = ax[0, 0].hist2d(row_plot1.compressed().ravel(),
                             d1[yvalue].compressed().ravel(),
                             **kw0)
        # kw0.update({"range": [[0, 4096] , [0.9, 0.9 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [0.7, 0.7 + 3.*ywide]],})
        p2 = ax[0, 1].hist2d(row_plot2.compressed().ravel(),
                             d2[yvalue].compressed().ravel(),
                             **kw0)
        #
        # kw0.update({"range": [[0, 4096] , [0.84, 0.84 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [1.02, 1.02 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [0.88, 0.88 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [0.4, 0.4 + 3.*ywide]],})
        p3 = ax[0, 2].hist2d(row_plot3.compressed().ravel(),
                             d3[yvalue].compressed().ravel(),
                             **kw0)
        #
        # kw0.update({"range": [[0, 4096] , [0.72, 0.72 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [1.06, 1.06 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [0.79, 0.79 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [-0.1, -0.1 + 5.*ywide]],})
        p4 = ax[0, 3].hist2d(row_plot4.compressed().ravel(),
                             d4[yvalue].compressed().ravel(),
                             **kw0)
        #
        # kw0.update({"range": [[0, 4096] , [0.45, 0.45 + 3.*ywide]],})
        # kw0.update({"range": [[0, 4096] , [1.09, 1.09 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [-1.8, -1.8 + 15.*ywide]],})
        p5 = ax[0, 4].hist2d(row_plot5.compressed().ravel(),
                             d5[yvalue].compressed().ravel(),
                             **kw0)
        # Normal processing'
        # kw0.update({"range": [[0, 4096] , [0.92, 0.92 + ywide]],})
        # kw0.update({"range": [[0, 4096] , [0.94, 0.94 + ywide]],})
        kw0.update({"range": [[0, 4096] , [0.98, 0.98 + ywide]],})
        p6 = ax[0, 5].hist2d(row_plot6.compressed().ravel(),
                             d6[yvalue].compressed().ravel(),
                             **kw0)
        # Some statistics histogram
        # dtype=[('norm', '<f4'), ('med', '<f4'), ('avg', '<f4'),
        # ('med_n', '<f4'), ('avg_n', '<f4'), ('rms_n', '<f4'),
        # ('unc_n', '<f4'), ('mad_n', '<f4')])
        parameter = "med_n"
        xnorm1 = d1[parameter].compressed().ravel() 
        xnorm2 = d2[parameter].compressed().ravel()
        xnorm3 = d3[parameter].compressed().ravel()
        xnorm4 = d4[parameter].compressed().ravel()
        xnorm5 = d5[parameter].compressed().ravel()
        xnorm6 = d6[parameter].compressed().ravel()
        kw_hist = {
            "bins": "auto",#""doane",# "auto",
            "histtype": "stepfilled",
            "align": "mid",
            "orientation": "vertical",
            "log": False,
            "color": "gray",
            }
        h1 = ax[1, 0].hist(xnorm1, **kw_hist)
        h2 = ax[1, 1].hist(xnorm2, **kw_hist)
        h3 = ax[1, 2].hist(xnorm3, **kw_hist)
        h4 = ax[1, 3].hist(xnorm4, **kw_hist)
        h5 = ax[1, 4].hist(xnorm5, **kw_hist)
        h6 = ax[1, 5].hist(xnorm6, **kw_hist)
        colors = ["gold", "lime", "cyan", "red", "gold"]
        for ind, qx in enumerate([5, 25, 50, 75, 95]):
            ax[1, 0].axvline(x=np.percentile(xnorm1, qx), c=colors[ind],
                            label="{0}%".format(qx))
            ax[1, 1].axvline(x=np.percentile(xnorm2, qx), c=colors[ind],
                            label="{0}%".format(qx))
            ax[1, 2].axvline(x=np.percentile(xnorm3, qx), c=colors[ind],
                            label="{0}%".format(qx))
            ax[1, 3].axvline(x=np.percentile(xnorm4, qx), c=colors[ind],
                            label="{0}%".format(qx))
            ax[1, 4].axvline(x=np.percentile(xnorm5, qx), c=colors[ind],
                            label="{0}%".format(qx))
            ax[1, 5].axvline(x=np.percentile(xnorm6, qx), c=colors[ind],
                            label="{0}%".format(qx))
        # Setup plot limits for histograms
        # ax[1, 0].set_xlim([60, 200])
        # ax[1, 1].set_xlim([60, 200])
        # ax[1, 2].set_xlim([60, 200])
        # ax[1, 3].set_xlim([60, 200])
        # ax[1, 4].set_xlim([-100, 100])
        ax[1, 0].set_xlim([np.percentile(xnorm1, 2), np.percentile(xnorm1, 98)])
        ax[1, 1].set_xlim([np.percentile(xnorm2, 2), np.percentile(xnorm2, 98)])
        ax[1, 2].set_xlim([np.percentile(xnorm3, 2), np.percentile(xnorm3, 98)])
        ax[1, 3].set_xlim([np.percentile(xnorm4, 2), np.percentile(xnorm4, 98)])
        ax[1, 4].set_xlim([np.percentile(xnorm5, 2), np.percentile(xnorm5, 98)])
        ax[1, 5].set_xlim([np.percentile(xnorm6, 2), np.percentile(xnorm6, 98)])
        #
        # ax[1, 0].set_xlim([70, 220])
        # ax[1, 1].set_xlim([70, 220])
        # ax[1, 2].set_xlim([70, 220])
        # Setup legend for histograms
        handles_ax1, labels_ax1 = ax[1, 0].get_legend_handles_labels()
        handles_ax2, labels_ax2 = ax[1, 1].get_legend_handles_labels()
        handles_ax3, labels_ax3 = ax[1, 2].get_legend_handles_labels()
        handles_ax4, labels_ax4 = ax[1, 3].get_legend_handles_labels()
        handles_ax5, labels_ax5 = ax[1, 4].get_legend_handles_labels()
        handles_ax6, labels_ax6 = ax[1, 5].get_legend_handles_labels()
        kw_lab = {
            "loc": "upper right",
            "ncol": 1,
            "markerscale": 2.1,
            "frameon": True,
            "fancybox": True,
            "framealpha": 0.5,
            "fontsize": 9,}
        ax[1, 0].legend(handles_ax1, labels_ax1, **kw_lab)
        ax[1, 1].legend(handles_ax2, labels_ax2, **kw_lab)
        ax[1, 2].legend(handles_ax3, labels_ax3, **kw_lab)
        ax[1, 3].legend(handles_ax4, labels_ax4, **kw_lab)
        ax[1, 4].legend(handles_ax5, labels_ax5, **kw_lab)
        ax[1, 5].legend(handles_ax6, labels_ax6, **kw_lab)
        # Colorbars
        # Set an UNIQUE range of values, 
        if False:
            kw_cbar = {"orientation": "horizontal", "aspect": 50}
            cb1 = plt.colorbar(p1[-1], ax=ax[0, 0], **kw_cbar)
            cb1.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
            cb1.set_label("N")
            cb2 = plt.colorbar(p2[-1], ax=ax[0, 1], **kw_cbar)
            cb2.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
            cb2.set_label("N")
            cb3 = plt.colorbar(p3[-1], ax=ax[0, 2], **kw_cbar)
            cb3.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
            cb3.set_label("N")
            cb4 = plt.colorbar(p4[-1], ax=ax[0, 3], **kw_cbar)
            cb4.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
            cb4.set_label("N")
            cb5 = plt.colorbar(p5[-1], ax=ax[0, 4], **kw_cbar)
            cb5.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
            cb5.set_label("N")
            cb6 = plt.colorbar(p6[-1], ax=ax[0, 5], **kw_cbar)
            cb6.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
            cb6.set_label("N")
            # Fontsize of the numbers!
            cb1.ax.tick_params(labelsize=7)
            cb2.ax.tick_params(labelsize=7)
            cb3.ax.tick_params(labelsize=7)
            cb4.ax.tick_params(labelsize=7)
            cb5.ax.tick_params(labelsize=7)
            cb6.ax.tick_params(labelsize=7)
        else:
            kw_cbar = {"orientation": "horizontal", "aspect": 20}
            # cb6 = plt.colorbar(p6[-1], ax=ax[0, 5], **kw_cbar)
            # cb6.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
            # cb6.set_label("N")
            # Fontsize of the numbers!
            # cb6.ax.tick_params(labelsize=7)
            cb6 = plt.colorbar(p6[-1], ax=ax[0, 5], **kw_cbar)
            cb6.patch.set_facecolor((0.2, 0.2, 0.2, 1.0))
            cb6.set_label("N")
            cb6.ax.tick_params(labelsize=7)
        # Setup ticks
        majorLocator_x0 = MultipleLocator(1024)
        majorFormatter_x0 = FormatStrFormatter("%d")
        minorLocator_x0 = MultipleLocator(128)
        ax[0, 0].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 0].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 0].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 1].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 1].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 1].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 2].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 2].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 2].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 3].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 3].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 3].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 4].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 4].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 4].xaxis.set_major_formatter(majorFormatter_x0)
        ax[0, 5].xaxis.set_minor_locator(minorLocator_x0)
        ax[0, 5].xaxis.set_major_locator(majorLocator_x0)
        ax[0, 5].xaxis.set_major_formatter(majorFormatter_x0)
        # For CCD3 plt.suptitle("CCD {0}, issue Amp".format(CCD))
        ttext = "CCD {0}, amp {1}.".format(CCD, amp)
        ttext += " Comparison of xtalk --version 2-- coeffs variations."
        ttext += " NO OVERSCAN SUBTRACTED FOR CCD02."
        ttext += "\nMedian normalized imgs. (victim-source)"
        if not (region is None):
            ttext += "   >>> REGION: {0}".format(region)
        if not (column is None):
            ttext += "   >>> COLUMNS: {0}".format(",".join(map(str, column)))
        plt.suptitle(ttext)
        # Labels
        kw_lab = {"fontsize": 10, "color": "blue",}
        # ax[0, 0].set_title("vers2. (3B-2B)/20", **kw_lab)
        # ax[0, 1].set_title("vers2.A (3B-2B)*0", **kw_lab)
        # ax[0, 2].set_title("vers2.B (3B-2B)*20", **kw_lab)
        # ax[0, 3].set_title("vers2.C (3B-2B)*40", **kw_lab)
        # ax[0, 4].set_title("vers2.D (3B-2B)*100", **kw_lab)
        # ax[0, 5].set_title("normal process", **kw_lab)
        # ax[0, 0].set_title("vers2.E (3B-2B)*1.7", **kw_lab)
        # ax[0, 1].set_title("vers2.F (3B-2B)*2.4", **kw_lab)
        # ax[0, 2].set_title("vers2.G (3B-2B)*5", **kw_lab)
        # ax[0, 3].set_title("vers2.H (3B-2B)*10", **kw_lab)
        # ax[0, 4].set_title("vers2.B (3B-2B)*20. 3X HEIGHT", **kw_lab) 
        # ax[0, 0].set_title("vers2.G (3B-2B)*5. 3X HEIGHT", **kw_lab)
        # ax[0, 1].set_title("vers2.H (3B-2B)*10. 3X HEIGHT", **kw_lab)
        # ax[0, 2].set_title("vers2.B (3B-2B)*20. 3X HEIGHT", **kw_lab)
        # ax[0, 3].set_title("vers2.C (3B-2B)*40. 5X HEIGHT", **kw_lab)
        # ax[0, 4].set_title("vers2.D (3B-2B)*100. 15X HEIGHT", **kw_lab)
        ax[0, 0].set_title("vers2.G (3B-2B)*5", **kw_lab)
        ax[0, 1].set_title("vers2.H (3B-2B)*10", **kw_lab)
        ax[0, 2].set_title("vers2.B (3B-2B)*20", **kw_lab)
        ax[0, 3].set_title("vers2.C (3B-2B)*40", **kw_lab)
        ax[0, 4].set_title("vers2.D (3B-2B)*100", **kw_lab)
        ax[0, 5].set_title("unaltered xtalk", **kw_lab)
        kw_lab.update({"color": "red"})
        ax[0, 0].set_ylabel(r"$med({box})$")
        ax[0 ,0].set_xlabel("rows")
        ax[0, 1].set_xlabel("rows")
        ax[0, 2].set_xlabel("rows")
        ax[0, 3].set_xlabel("rows")
        ax[0, 4].set_xlabel("rows")
        ax[0, 5].set_xlabel("rows")
        ax[1, 0].set_ylabel("N")
        ax[1, 0].set_xlabel(r"$med({box})$")
        ax[1, 1].set_xlabel(r"$med({box})$")
        ax[1, 2].set_xlabel(r"$med({box})$")
        ax[1, 3].set_xlabel(r"$med({box})$")
        ax[1, 4].set_xlabel(r"$med({box})$")
        ax[1, 5].set_xlabel(r"$med({box})$")
        plt.subplots_adjust(left=0.05, bottom=0.06, right=0.98, top=0.90,
                            wspace=0.26, hspace=0.35)
        if suffix is None:
            suffix = "pid" + str(os.getpid())
        if saveplot:
            if (region is None):
                outname = "profileComp_vers2_{0:02}_{1}.png".format(ccd,suffix)
            else:
                outname = "profileComp_vers2_{2}_{0:02}_{1}.png".format(ccd,
                                                                        suffix,
                                                                        region)
            plt.savefig(outname, dpi=400, facecolor="w", edgecolor="w",
                    orientation="portrait", papertype=None, format="png",
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    frameon=None)
        plt.show()


class Explore():
    """ Class to harbor methods for explore sub-populations in the row profile
    """
    def __init__(self, ccd=None, amp=None, saveplot=True,
                 suffix=None, region=None):
        self.ccd = ccd
        self.amp = amp
        self.saveplot = saveplot
        self.suffix = suffix
        self.region = region

    def expl1(self):
        """ Method to explore different sub-populations from the row profile
        distribution
        """
        #
        # Get the data
        #
        if (self.region is None):
            fnm1 = "stat_{0:02}_16x128_medN_v2.npy".format(CCD)
            fnm2 = "stat_{0:02}_16x128_medN_v2.A.npy".format(CCD)
            fnm3 = "stat_{0:02}_16x128_medN_v2.B.npy".format(CCD)
            fnm4 = "stat_{0:02}_16x128_medN_v2.C.npy".format(CCD)
            fnm5 = "stat_{0:02}_16x128_medN_v2.D.npy".format(CCD)
            fnm6 = "stat_{0:02}_16x128_medN_normal.npy".format(CCD)
        else:
            fnm1 = "stat_{0:02}_16x128_medN_v2_{1}.npy".format(CCD, region)
            fnm2 = "stat_{0:02}_16x128_medN_v2.A_{1}.npy".format(CCD, region)
            fnm3 = "stat_{0:02}_16x128_medN_v2.B_{1}.npy".format(CCD, region)
            fnm4 = "stat_{0:02}_16x128_medN_v2.C_{1}.npy".format(CCD, region)
            fnm5 = "stat_{0:02}_16x128_medN_v2.D_{1}.npy".format(CCD, region)
            fnm6 = "stat_{0:02}_16x128_medN_normal.npy".format(CCD)
        data = [fnm1, fnm2, fnm3, fnm4, fnm5, fnm6]
        # Retrieve only one amplifier, masked arrays
        tmp1, tmp2, tmp3 = np.load(fnm1), np.load(fnm2), np.load(fnm3)
        tmp4, tmp5, tmp6 = np.load(fnm4), np.load(fnm5), np.load(fnm6)
        row_plot = None
        d1, row_plot1 = Compare().sel_amplifier(tmp1, amp=self.amp, 
                                                ccdnum=self.ccd)
        d2, row_plot2 = Compare().sel_amplifier(tmp2, amp=self.amp, 
                                                ccdnum=self.ccd)
        d3, row_plot3 = Compare().sel_amplifier(tmp3, amp=self.amp, 
                                                ccdnum=self.ccd)
        d4, row_plot4 = Compare().sel_amplifier(tmp4, amp=self.amp, 
                                                ccdnum=self.ccd)
        d5, row_plot5 = Compare().sel_amplifier(tmp5, amp=self.amp, 
                                                ccdnum=self.ccd)
        d6, row_plot6 = Compare().sel_amplifier(tmp6, amp=self.amp, 
                                                ccdnum=self.ccd)
        #
        # Get the additional info
        #
        if (self.region is None):
            i1 = "info_{0:02}_16x128_medN_v2.npy".format(CCD)
            i2 = "info_{0:02}_16x128_medN_v2.A.npy".format(CCD)
            i3 = "info_{0:02}_16x128_medN_v2.B.npy".format(CCD)
            i4 = "info_{0:02}_16x128_medN_v2.C.npy".format(CCD)
            i5 = "info_{0:02}_16x128_medN_v2.D.npy".format(CCD)
            i6 = "info_{0:02}_16x128_medN_normal.npy".format(CCD)
        else:
            i1 = "info_{0:02}_16x128_medN_v2_{1}.npy".format(CCD, region)
            i2 = "info_{0:02}_16x128_medN_v2.A_{1}.npy".format(CCD, region)
            i3 = "info_{0:02}_16x128_medN_v2.B_{1}.npy".format(CCD, region)
            i4 = "info_{0:02}_16x128_medN_v2.C_{1}.npy".format(CCD, region)
            i5 = "info_{0:02}_16x128_medN_v2.D_{1}.npy".format(CCD, region)
            i6 = "info_{0:02}_16x128_medN_normal.npy".format(CCD)
        info = [i1, i2, i3, i4, i5, i6]
        # Retrieve only one amplifier, masked arrays
        t1, t2, t3 = np.load(i1), np.load(i2), np.load(i3)
        t4, t5, t6 = np.load(i4), np.load(i5), np.load(i6)
        #
        # Do the math
        #
        # For each exposure, look at the cummulative sum, and associate it
        # with its observational data
        x0, x1, x2, x3 = [], [], [], []
        for idx, k in enumerate(d1):
            x0.append(t1[idx]["nite"])
            aux = d1["med_n"][idx].compressed().ravel()
            x1.append(np.sum(aux))
            x2.append(np.median(aux))
            x3.append(np.median(np.abs(aux - np.median(aux))))
        df1 = pd.DataFrame({"nite":x0, "sum":x1, "median": x2, "mad": x3})
        #
        y0, y1, y2 = [], [], []
        for n in np.unique(df1["nite"].values[:]):
            y0.append(n[0])
            y1.append(np.median(df1.loc[df1["nite"]==n[0], "sum"].values))
            y2.append(np.median(df1.loc[df1["nite"]==n[0], "median"].values))
        #
        plt.close("all")
        fig, ax = plt.subplots(2, 2, figsize=(8, 8))
        #
        # Sort by expnumand transform nite to a plottable quantity
        xdate = [datetime.datetime.strptime(str(date[0]), "%Y%m%d")
                for date in np.sort(df1["nite"])]
        xdate2 = [datetime.datetime.strptime(str(date), "%Y%m%d")
                  for date in y0]
        #
        ax[0, 0].scatter(xdate, df1["sum"], c="deepskyblue", marker=".",
                         s=20)
        ax[0, 1].scatter(xdate, df1["sum"], c="deepskyblue", marker=".",
                         s=20)
        ax[1, 0].scatter(xdate, df1["median"], c="goldenrod", marker="o", s=20)
        ax[1, 1].scatter(xdate, df1["median"], c="goldenrod", marker="o", s=20)
        ax[0, 1].scatter(xdate2, y1, c="k", marker="^", s=10)
        ax[1, 1].scatter(xdate2, y2, c="k", marker="^", s=10)
        # from 0.95 to 1.05
        ax[0, 1].set_ylim([0.97*1792, 1.01*1792])
        ax[1, 1].set_ylim([0.97, 1.01])
        # 
        months = MonthLocator()
        days = DayLocator(interval=10)
        dateFmt = DateFormatter("%b/%d/%y")
        months2 = MonthLocator()
        days2 = DayLocator()
        dateFmt2 = DateFormatter("%b/%d/%y")
        #
        for xx in [ax[0,0], ax[0,1], ax[1,0], ax[1,1]]:
            xx.xaxis.set_major_locator(months)
            xx.xaxis.set_major_formatter(dateFmt)
            # xx.xaxis.set_minor_locator(days)
            xx.autoscale_view()
        kw_grid1 = {
            "color": "lightgray",
            "linestyle": "dotted",
            "dash_capstyle": "round",
            "alpha": 0.7,}
        for XX in [ax[0,0], ax[0,1], ax[1,0], ax[1,1]]:
            XX.grid(**kw_grid1)
            xlabels_ax1 = XX.get_xticklabels()
            plt.setp(xlabels_ax1, rotation=30, fontsize=10)
        # For the titles
        ax[0, 0].set_title("Sum of boxes (normalized) median, N values")
        ax[0, 1].set_title("Zoom in")
        ax[1, 0].set_title("Median of the boxes, 1 value per exposure")
        ax[1, 1].set_title("Zoom in")
        #
        plt.subplots_adjust(left=0.12, bottom=0.09, right=0.90, top=0.95,
                        wspace=0.20, hspace=0.35)
        if self.suffix is None:
            self.suffix = "pid" + str(os.getpid())
        if self.saveplot:
            if (self.region is None):
                outname = "profile_explore_{0:02}_{1}.png".format(self.ccd, 
                                                          self.suffix)
            else:
                outname = "profile_explore_{2}_{0:02}_{1}.png".format(self.ccd,
                                                              self.suffix,
                                                              self.region)
            plt.savefig(outname, dpi=400, facecolor="w", edgecolor="w",
                    orientation="portrait", papertype=None, format="png",
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    frameon=None)
        plt.show()
        return True


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
    #
    # Amplifier: use A and B instead of 1 and 0, because it gave me some 
    # headaches
    
    CCD = 3
    amplifier = "A"

    C = Compare()
    C.side4(ccd=CCD, amp=amplifier, saveplot=False, 
            region="reg3",
            suffix="xtmodif_medN_amp{0}_col0_6_opt2".format(amplifier),
            column=[0, 6])
    
    exit()
    
    for iter_col in range(8):
        C.side4(ccd=CCD, amp=amplifier, saveplot=True, 
                region="reg2",
                suffix="xtmodif_medN_amp{1}_col{0}".format(iter_col, amplifier),
                column=[iter_col])
    

    E = Explore(ccd=CCD, amp=amplifier, suffix=str(os.getpid()), 
                saveplot=True,
                region=None)
    E.expl1()

    
    #C.side3(ccd=CCD, amp=amplifier, saveplot=False, suffix=str(os.getpid()))
    #
    # C.side2(ccd=CCD, amp=0, saveplot=True, suffix="AmpB_set2")
    # C.side2(ccd=CCD, amp=1, saveplot=True, suffix="AmpA_set2")
    # C.side1(ccd=CCD, amp=1, saveplot=True, suffix="AmpA_set2")
    # C.side1(ccd=CCD, amp=0, saveplot=True, suffix="AmpB_set2")
    # C.side2(ccd=CCD, amp=1, saveplot=False, suffix="AmpA_set1")
    # C.side1(ccd=CCD, amp=0, saveplot=False, suffix="AmpB_set1")


    # Histo = Tide(stat="stat_{0:02}_16x128_y4e1_medN.npy".format(CCD),
    #         db="info_{0:02}_16x128_y4e1_medN.npy".format(CCD))
    # Histo = Tide(stat="stat_{0:02}_xtalked_16x128_medN.npy".format(CCD),
    #         db="info_{0:02}_xtalked_16x128_medN.npy".format(CCD))
    # Histo = Tide(stat="stat_{0:02}_noOsc_16x128_medN.npy".format(CCD),
    #          db="info_{0:02}_noOsc_16x128_medN.npy".format(CCD))
    # Histo.histo_param(ccd=CCD, suffix="noOsc_16x128_medN")

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
