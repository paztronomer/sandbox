""" Script to evaluate the circular symmetry of sources 
"""

import os
import glob
import time
import argparse
import logging
import numpy as np
from scipy import interpolate
from scipy import stats
from skimage.measure import compare_ssim
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize, ZScaleInterval)
import fitsio

def draw_ellipse():
    # This method was saved for simplicity but is not working
    # Plot the stamp using min/max coordinates, adding the ellipse
    fig, ax = plt.subplots(1, 2, figsize=(4, 2))
    im = ax[0].imshow(stamp, origin='lower', cmap='gray_r', 
                      interpolation='none')
    # Plot centroid
    ax[0].plot(xcnt, ycnt, 'ro')
    # Plot ellipse
    ellipse = mpatches.Ellipse([xcnt, ycnt], 
                               2 * a, 
                               2 * b, 
                               angle=theta,
                               edgecolor='yellow', 
                               facecolor=None, 
                               fill=False, 
                               linewidth=1.5)
    circle = mpatches.Circle([xcnt, ycnt], 
                              radius=2 * a, 
                              edgecolor='tomato',
                              linestyle='-',
                              facecolor=None, 
                              fill=False, 
                              linewidth=1.5)
    ax[0].add_artist(ellipse)
    ax[0].add_artist(circle)
    # Draw lines for circle section
    ax[0].plot([xcnt, xl], [ycnt, yl_up], '-', lw=1, color='lime')
    ax[0].plot([xcnt, xl], [ycnt, yl_dw], '-', lw=1, color='lime')
    ax[0].plot([xcnt, xr], [ycnt, yr_up], '-', lw=1, color='lime')
    ax[0].plot([xcnt, xr], [ycnt, yr_dw], '-', lw=1, color='lime')
    # Set info
    txt0 = '({0}, {1}), a/b={2:.2f}'.format(xcnt, ycnt, a / b)
    ax[0].set_title(txt0)
    plt.show()

def open_fits(fnm1, fnm2, ext1=1, ext2=2):
    """ Open both files, the image file and the catalog
    """
    with fitsio.FITS(fnm1) as f1:
        x1 = f1[ext1].read()
        h1 = f1[ext1].read_header()
    with fitsio.FITS(fnm2) as f2:
        x2 = f2[ext2].read()
        h2 = f2[ext2].read_header()
    return (x1, h1), (x2, h2)

def table_catalog(cat, h_cat):
    """ Receives catalog array, and header from the fits. It returns a pandas
    dataframe with the lowercase columns. It contains all the fields and the 
    expanded APER fields
    Inputs
    - cat: catalog fits file
    - h_cat: header for catalog fits file 
    """
    # Construct the catalog
    # Each set of entries is divided in 4: TTYPEnn, TFORMnn, TUNITnn, TDISPnn
    # Warning: APER quantities has an array stored intead of a single value
    aux_colnames = []
    aux_data = []
    aux_keys= [h_cat[i].strip().lower() 
               for i in h_cat.keys() if 'TTYPE' in i]
    # We need more column names than the actual number, because APER keys
    # need to be expanded. We will use the first element of the structured
    # in combination with the list of keywords previously defined
    for idx_k, k in enumerate(cat[0]): 
        # If the element is an array, expand its name
        if (isinstance(k, np.ndarray)):
            for idx_aper, aper in enumerate(k):
                tmp_nm = '{0}_{1}'.format(aux_keys[idx_k], 
                                          idx_aper + 1)
                aux_colnames.append(tmp_nm)
        else:
            aux_colnames.append(aux_keys[idx_k])
    # Use all rows to get expanded rows. The overall datatype is ndarray, but 
    # the rows dtype is void. This dtype (void) can harbor different dytpes
    for row in cat:
        aux_row = []
        for k in row:
            if (isinstance(k, np.ndarray)):
                aux_row += list(k)
            else:
                aux_row.append(k)
        aux_data.append(aux_row)
    # Construct the dataframe
    df_cat = pd.DataFrame(aux_data, columns=aux_colnames)
    return df_cat

def region_lr(cnt_coo, a, 
              rad_factor=1,
              angle_lr=[np.pi, 0], 
              alpha=10, 
              sampling_factor=3,
              img=None,
              do_plot=False,):
    """ Method to get data from 2 circular regions at the left and right
    or a centroid. Returns 2 masked arrays, one for the left and other for
    the right side of the centroid.
    The values of the fine grid images have its values normalized to sum the
    initial value of its parent pixel.
    This method works well if:
    * borders of circular regions are not +/- pi/2
    * both regions can be defined with only 2 straight lines
    
    To refine: do not use the stamp but do calculations directly on the CCD
    section
    
    Inputs
    - cnt_coo: centroid, where coordinates are [x_image, y_image]
    - radius: adjusted radio of the object
    - rad_factor: factor to be multiplied for the major axis value of the
    elliptical adjust to the source
    - angle_lr: left, right reference angles to create the circular regions 
    - alpha: angle in degrees to set a range above and below the left/right
    reference angles
    - sampling_factor: factor to supersample in both direction
    - img: image to be displayed for visualization
    - do_plot: boolean to set the plotting of the stamp
    """
    #
    # Creation of the stamp: interpolation and selected regions
    #
    # Get the circle sectors, for left and right of the centroid
    # Radius for statistics
    radius = rad_factor * a
    # Circular section 
    # Upper and lower angles
    theta_l_up = angle_lr[0] - np.deg2rad(alpha)
    theta_l_dw = angle_lr[0] + np.deg2rad(alpha)
    theta_r_up = angle_lr[1] + np.deg2rad(alpha)
    theta_r_dw = angle_lr[1] - np.deg2rad(alpha)
    # (x, y) coordinates of limit positions, from polar coordinates
    # These coordinates are centered at the origin
    xl = radius * np.cos(theta_l_up)
    yl_up, yl_dw = radius * np.sin(theta_l_up), radius * np.sin(theta_l_dw)
    xr = radius * np.cos(theta_r_up)
    yr_up, yr_dw = radius * np.sin(theta_r_up), radius * np.sin(theta_r_dw)
    #
    # Stamp centroid and borders. The stamps is only used as auxiliary for
    # selecting circular regions
    # Note yr_dw == yl_dw and yr_up == yl_up for symmetric sections
    xc_s, yc_s = radius, radius
    xl_s, xr_s = 0, 2. * radius
    yl_dw_s, yr_dw_s = 0, 0 
    yl_up_s, yr_up_s = 2. * radius, 2. * radius
    # xc_s, yc_s = np.abs(xl), np.abs(yl_dw)
    # xl_s, yl_dw_s, yl_up_s = 0, 0, np.ptp([yl_dw, yl_up])
    # xr_s, yr_dw_s, yr_up_s = np.ptp([xl, xr]), 0, np.ptp([yr_dw, yr_up])
    if not True:
        print(xc_s, yc_s, xl_s, yl_dw_s, yl_up_s, xr_s, yr_dw_s, yr_up_s)
    #
    # CCD region centroid and borders, use radius to define region
    # C=Object centroid is at the image center
    xc_ccd, yc_ccd = cnt_coo
    xl_ccd, xr_ccd = xc_ccd - radius , xc_ccd + radius
    yl_dw_ccd, yl_up_ccd = yc_ccd - radius, yc_ccd + radius
    yr_dw_ccd, yr_up_ccd = yc_ccd - radius, yc_ccd + radius
    # xl_ccd, xr_ccd = xl + xc_ccd, xr + xc_ccd
    # yl_dw_ccd, yl_up_ccd = yl_dw + yc_ccd, yl_up + yc_ccd
    # yr_dw_ccd, yr_up_ccd = yr_dw + yc_ccd, yr_up + yc_ccd
    # The precission value 0.001 was defined arbitrary
    if (((yl_dw - yr_dw) < 0.001) and ((yl_up - yr_up) < 0.001)):
        pass
    else:
        print('ERROR')
        logging.error('Circular sections are not symmetric')
        exit(1)
    #
    # Circular section auxiliary functions
    #
    # Get the pixels inside the circle sections: 2 lines needs to be adjusted
    # Outputs from linregress are: m, b, r_value, p_value, std_err
    # Then we'll evaluate each one of the coordinates to see if they belong 
    # to circle left/right sections
    # These 2 lines are fitted on the stamp coordinate system
    # f = interpolate.interp1d([xl_s, xc_s, xr_s], 
    #                          [yl_up_s, yc_s, yr_dw_s], 
    #                          kind='linear',
    #                          bounds_error=False,
    #                          fill_value='extrapolate')
    # g = interpolate.interp1d([xl_s, xc_s, xr_s], 
    #                          [yl_dw_s, yc_s, yr_up_s], 
    #                          kind='linear',
    #                          bounds_error=False,
    #                          fill_value='extrapolate')
    # Linear regression
    m1, b1, r1, p1, stderr1 = stats.linregress([xl_s, xc_s, xr_s], 
                                               [yl_up_s, yc_s, yr_dw_s],)
    m2, b2, r2, p2, stderr2 = stats.linregress([xl_s, xc_s, xr_s], 
                                               [yl_dw_s, yc_s, yr_up_s],)
    f = lambda x: m1 * x + b1
    g = lambda x: m2 * x + b2
    #
    # Create a grid with higher resolution
    #
    # Increase the number of points of the grid to each pixel be NxN as defined
    # by sampling_factor
    # Use left side y-coordinate to define the grid
    xmin_grid, xmax_grid = int(np.floor(xl_ccd)), int(np.ceil(xr_ccd))
    ymin_grid, ymax_grid = int(np.floor(yl_dw_ccd)), int(np.ceil(yl_up_ccd))
    # add some range to the y-axis
    # ymin_grid -= 5
    # ymax_grid += 5
    # Number of total boxes the higher resolution wil contain
    Nx = (sampling_factor) * np.ptp([xmax_grid, xmin_grid]) + 1
    Ny = (sampling_factor) * np.ptp([ymax_grid, ymin_grid]) + 1
    # Define the supersampled grid
    yy, xx = np.meshgrid(np.linspace(ymin_grid, ymax_grid, Ny), 
                         np.linspace(xmin_grid, xmax_grid, Nx),
                         sparse=False, indexing='ij')
    #
    # Interpolate grid values
    #
    # 1) Define initial grid from the image, and transform it to set of 
    # coordinates
    yy_ini, xx_ini = np.mgrid[ymin_grid:ymax_grid, xmin_grid:xmax_grid]
    points = np.vstack([xx_ini.ravel(), yy_ini.ravel()])
    values = img[ymin_grid:ymax_grid, xmin_grid:xmax_grid].flatten()
    # 2) Using the positions and values from the original image, flesh out the
    # finer grid. Nearest value will be assumed
    interp_img = interpolate.griddata(
        points.T, 
        values, 
        (xx, yy), 
        method='nearest'
    )
    # 3) Normalize by the number of subpixels each pixel is divided, then I'll
    # not create additional values
    interp_img /= np.power(sampling_factor, 2)
    #
    # NOTE: xx, xx_ini, yy, yy_ini have the same range of values, based on the 
    # original CCD positions. The interpolated grid (interp_img) has positions
    # shifted to start at the origin.
    if False:
        plt.imshow(interp_img)
        plt.show()
    
    # To get the region of interest from the original image
    # img[ymin_grid:ymax_grid, xmin_grid:xmax_grid]

    #
    # Masks for the circular regions
    #
    # Select the points from the fine grid belonging to the region of interest 
    # and create a mask. Make sure to mask using the correct set of coordinates
    #
    # Circle function
    h = lambda x, y, X, Y: np.sqrt(np.power(x - X, 2.) + np.power(y - Y, 2.)) 
    # Define set of coordinates for the stamp. Lowe corner is the origin
    # xx_s, yy_s = xx - xl_ccd, yy - yl_dw_ccd
    xx_s, yy_s = xx - xmin_grid, yy - ymin_grid
    # 1) Left/right half circle: use the set of coordinates from the
    # grid interpolation process to define the pixels positions
    # Get the mask
    zz_s = h(xx_s, yy_s, xc_s, yc_s) 
    interp_circle_msk = zz_s > radius
    #
    if not True:
        plt.plot(xx_s, yy_s, 'g.')
        plt.show()
    # Apply selection to left/right halves. Need to combine masks
    aux_c01 = np.ma.masked_where(interp_circle_msk, xx_s) 
    c01 = aux_c01 > xc_s
    aux_c02 = np.ma.masked_where(interp_circle_msk, xx_s) 
    c02 = aux_c02 < xc_s
    circ_l = np.ma.masked_where(c01, interp_img)
    circ_r = np.ma.masked_where(c02, interp_img)
    circ_l_msk = np.ma.getmask(circ_l)
    circ_r_msk = np.ma.getmask(circ_r)
    #
    if not True:
        plt.imshow(circ_l)
        plt.show()
    """
    # Circle function
    h = lambda x, y: np.sqrt(np.power(x, 2.) + np.power(y, 2.)) 
    # Define set of coordinates for the stamp
    xx_s, yy_s = xx - xl_ccd, yy - yl_dw_ccd
    # Move to the center of the stamp to have coordinates centered at zero
    xx_s -= xc_s
    yy_s -= yc_s
    # 1) Left/right half circle: use the set of coordinates from the
    # grid interpolation process to define the pixels positions
    # Get the mask
    zz_s = h(xx_s, yy_s) #np.sqrt(np.power(xx_s, 2.) + np.power(yy_s, 2.))
    interp_circle_msk = zz_s > radius
    # Apply selection to left/right halves. Need to combine masks
    aux_c01 = np.ma.masked_where(interp_circle_msk, xx_s) 
    c01 = aux_c01 > 0
    aux_c02 = np.ma.masked_where(interp_circle_msk, xx_s) 
    c02 = aux_c02 < 0
    circ_l = np.ma.masked_where(c01, interp_img)
    circ_r = np.ma.masked_where(c02, interp_img)
    circ_l_msk = np.ma.getmask(circ_l)
    circ_r_msk = np.ma.getmask(circ_r)
    """
    #
    # 2) Left/right angular circle sections
    c03 = ~np.logical_and(yy_s <= f(xx_s), yy_s >= g(xx_s))
    c04 = ~np.logical_and(yy_s >= f(xx_s), yy_s <= g(xx_s))
    c05 = np.ma.masked_where(interp_circle_msk, c03)
    c06 = np.ma.masked_where(interp_circle_msk, c04)
    angle_l = np.ma.masked_where(c05, interp_img)
    angle_r = np.ma.masked_where(c06, interp_img)
    #
    if not True:
        plt.imshow(np.ma.masked_where(angle_l, interp_img))
        plt.imshow(np.ma.masked_where(angle_r, interp_img)) 
        plt.show()
    # Now: get stats from these regions!!!
    # Plotting the resampled grid for evaluation
    if do_plot:
        stamp = img[ymin_grid:ymax_grid, xmin_grid:xmax_grid]
        im_norm = ImageNormalize(stamp, 
                                 interval=ZScaleInterval(),
                                 stretch=SqrtStretch())
        im_norm2 = ImageNormalize(interp_img, 
                                  interval=ZScaleInterval(),
                                  stretch=SqrtStretch())
        fig, ax = plt.subplots(1, 2, figsize=(8, 4))
        kw = {'origin': 'lower', 'cmap': 'viridis',}
        # i = ax[0, 0].imshow(img, norm=im_norm, **kw)
        # ax[0, 0].plot(xcnt, ycnt, 'ro')
        i = ax[0].imshow(img, norm=im_norm, **kw)
        # ax[1, 0].scatter(xx, yy, marker='.', s=5, color='white', alpha=0.3)
        ax[0].plot(xcnt, ycnt, 'ro')
        # i2 = ax[0, 1].imshow(interp_img, norm=im_norm2, **kw)
        i2 = ax[1].imshow(interp_img, norm=im_norm2, **kw)
        # Some text for clarity
        ax[0].set_title('Original resolution')
        ax[1].set_title('Interpolated image')
        # Circumference and circular sections
        circle = mpatches.Circle([xcnt, ycnt], 
                                  radius=radius, 
                                  edgecolor='tomato',
                                  linestyle='-',
                                  facecolor=None, 
                                  fill=False, 
                                  linewidth=1.5)
        circle2 = mpatches.Circle([xcnt - xmin_grid, ycnt - ymin_grid], 
                                   radius=radius, 
                                   edgecolor='tomato',
                                   linestyle='-',
                                   facecolor=None, 
                                   fill=False, 
                                   linewidth=1.5)
        ax[0].add_artist(circle)
        ax[1].add_artist(circle2)
        # Draw lines for circle section
        ax[0].plot([xcnt, xl_ccd], [ycnt, yl_up_ccd], '-', lw=1, color='k')
        ax[0].plot([xcnt, xl_ccd], [ycnt, yl_dw_ccd], '-', lw=1, color='k')
        ax[0].plot([xcnt, xr_ccd], [ycnt, yr_up_ccd], '-', lw=1, color='k')
        ax[0].plot([xcnt, xr_ccd], [ycnt, yr_dw_ccd], '-', lw=1, color='k')
        ax[1].plot([xcnt - xmin_grid, xl_ccd - xmin_grid], 
                   [ycnt - ymin_grid, yl_up_ccd - ymin_grid], 
                   '-', lw=1, color='k')
        ax[1].plot([xcnt - xmin_grid, xl_ccd - xmin_grid], 
                   [ycnt - ymin_grid, yl_dw_ccd - ymin_grid], 
                   '-', lw=1, color='k')
        ax[1].plot([xcnt - xmin_grid, xr_ccd - xmin_grid], [ycnt - ymin, yr_up - ymin], 
                   '-', lw=1, color='k')
        ax[1].plot([xcnt - xmin, xr - xmin], [ycnt - ymin, yr_dw - ymin], 
                   '-', lw=1, color='k')
        # Limits
        ax[0].set_xlim([xmin, xmax])
        ax[0].set_ylim([ymin, ymax])
        # ax[1].set_xlim([xmin, xmax])
        # ax[1].set_ylim([ymin, ymax])
        plt.suptitle('Original vs interpolated. Plot works based on symmetry')
        plt.show()
    # Return only left/right
    return angle_l, angle_r, circ_l, circ_r

def masked_stat(ma_a):
    """ Receives a masked array and performs some basic statistics
    """
    ma_a = ma_a.compressed()
    mad = lambda x: np.median(np.abs(x - np.median(x)))
    # Add some statistical tests
    out = [np.mean(ma_a), np.median(ma_a), np.std(ma_a), 
           np.var(ma_a), mad(ma_a), np.ptp(ma_a), ma_a.size]
    return out

def aux_main(outname='2region_stat.csv',
             d0_range=[1, 4096], 
             d1_range=[1024, 2048 - 100],
             table=None,
             ccdlist=None,
             ftypes=['red_immask', 'cat_firstcut'],
             rpath='/archive_data/desarchive'):
    """ Run the processing
    Inputs
    - d0_range: select the entire range in long dimension
    - d1_range: select the right amplifier, cutting off a vertical band at the 
    border for the tapebumps
    - outname: output name for the CSV results table
    - table: filename of the CSV table containing t_eff, path, filetype
    - ccdlist: 
    """
    # ccd = 41
    # Call opening of files
    # x1 = "/archive_data/desarchive/OPS/firstcut/Y6N/20180914-r3563/"
    # x1 += "D00773758/p01/red/immask/D00773758_g_c41_r3563p01_immasked.fits.fz"
    # x2 = "/archive_data/desarchive/OPS/firstcut/Y6N/20180914-r3563/"
    # x2 += "D00773758/p01/red/immask/D00773758_g_c28_r3563p01_immasked.fits.fz" 
    # y1 = "/archive_data/desarchive/OPS/firstcut/Y6N/20180914-r3563/"
    # y1 += "D00773758/p01/cat/D00773758_g_c41_r3563p01_red-fullcat.fits"
    # y2 = "/archive_data/desarchive/OPS/firstcut/Y6N/20180914-r3563/"
    # y2 += "D00773758/p01/cat/D00773758_g_c28_r3563p01_red-fullcat.fits"
    # fnm1 = [x1, x2]
    # fnm2 = [y1, y2]
    # ccd = [ccd]

    t0 = time.time()
    # Open table, identify filetypes
    df = pd.read_csv(table)
    df.columns = df.columns.map(str.lower)
    # Sort by path
    df.sort_values(['path'], inplace=True)
    # Reset index
    df.reset_index(drop=True, inplace=True)
    
    # Create the file list for both filetypes
    fnm1, fnm2 = [], []
    if (len(ftypes) != 2):
        logging.error('More than 2 filetypes. Need to change a block')
        exit()
    for c in ccdlist:
        a = [glob.glob(os.path.join(rpath, x, '*_c{0:02}_*'.format(c)))
             for x in df.loc[df['filetype'] == ftypes[0], 'path']
        ]
        b = [glob.glob(os.path.join(rpath, x, '*_c{0:02}_*'.format(c)))
             for x in df.loc[df['filetype'] == ftypes[1], 'path']
        ]
        if ((len(a) > 0) and (len(b) == len(a))):
            fnm1 += a
            fnm2 += b
        else:
            logging.warning('CCD:{0} error locating files. Skip'.format(c))
    # Flatten the arrays, because glob.glob introduces an additional level of
    # nesting
    fnm1 = np.array(fnm1).flatten()
    fnm2 = np.array(fnm2).flatten()

    # Results list to construct the dataframe
    res = []

    for i in range(len(fnm1)):
        t_f0 = time.time()
        fits_img, fits_cat = open_fits(fnm1[i], fnm2[i])
        sci, h_sci = fits_img
        cat, h_cat = fits_cat
        # From the array and header get a dataframe of the catalog
        dfcat = table_catalog(cat, h_cat)
        # Remove unused catalog array and header
        del cat, h_cat
        # Select stars in the patch of interest, and with FLAGS=0
        c1 = (dfcat['x_image'] >= d1_range[0])
        c2 = (dfcat['x_image'] <= d1_range[1])
        c3 = (dfcat['y_image'] >= d0_range[0])
        c4 = (dfcat['y_image'] <= d0_range[1])
        c5 = (dfcat['flags'] == 0)
        cond1 = c1 & c2
        cond2 = c3 & c4 & c5
        # Dataframe of the centroids of interest
        # Note the centroids are decimal
        sel = dfcat.loc[cond1 & cond2]

        #
        # Experimentation
        #
        aux_val = dfcat['mag_aper_4'].values
        aux_val = aux_val[np.where(aux_val < 90)]
        q = np.percentile(aux_val, 10)
        # Add the magnitude cut
        sel = sel.loc[sel['mag_aper_4'] < q]
        #
        #
        #
        print(os.path.basename(fnm1[i]), len(dfcat.index), len(sel.index))

        # Iterate over all the objects in the region. For each of them 
        # calculate the left/right flux statistics
       
        cnt_x = sel['x_image'].values
        cnt_y = sel['y_image'].values
        a = sel['a_image'].values
        
        # Create lists to be saved into a dataframe
        expnum = h_sci['EXPNUM']
        ccdnum = h_sci['CCDNUM']
        nite = int(h_sci['NITE'])
        band = h_sci['BAND'].strip()
        region = '[{0}:{1},{2}:{3}]'.format(*d1_range, *d0_range)
        # Also save T_EFF 
        t_eff = df.loc[df['expnum'] == expnum, 't_eff'].values[0]

        for obj in range(cnt_x.size): 
            # Using the above parameters define 2 circular regions 
            results = region_lr([cnt_x[obj], cnt_y[obj]], a[obj], 
                                rad_factor=1,
                                alpha=10,
                                angle_lr=[np.pi, 0], 
                                img=sci,
                                do_plot=False,)
            # Uncompress masked arrays
            angle_l, angle_r, circ_l, circ_r = results
            #
            # Get the statistics for each region
            aux_angle_l = masked_stat(angle_l) + ['angle', 'L']
            aux_angle_r = masked_stat(angle_r) + ['angle', 'R']
            aux_circ_l = masked_stat(circ_l) + ['circ', 'L']
            aux_circ_r = masked_stat(circ_r) + ['circ', 'R']
            # Get fitted parameters for each object
            fpar = [cnt_x[obj], 
                    cnt_y[obj], 
                    sel['mag_aper_4'].values[obj],
                    sel['a_image'].values[obj], 
                    sel['b_image'].values[obj],
                    sel['theta_image'].values[obj], 
                    sel['number'].values[obj],]
            # Fill the output list
            tmp_data = [expnum, ccdnum, nite, band, t_eff, region]
            tmp_data += fpar
            for reg in [aux_angle_l, aux_angle_r, aux_circ_l, aux_circ_r]:
                t = tmp_data + reg
                res.append(t)
            
            # -----------------------------------------------------------------
            # Construct an quick assessment plot, based in basic statistics 
            if False:
                # Quick plot distribution of a, b
                fig, ax = plt.subplots(2,2)
                sel = dfcat.loc[dfcat['x_image'] > 1024]
                sel = sel.loc[sel['theta_image'].abs() < 10]
                sel0 = dfcat.loc[dfcat['x_image'] < 1024]
                sel0 = sel0.loc[sel0['theta_image'].abs() < 10]
                ax[0, 0].hist(sel['a_image'] / sel['b_image'], 
                              histtype='step', bins=20,
                              color='g')
                ax[0, 1].plot(sel['a_image'] / sel['b_image'], 
                              sel['theta_image'], 'go')
                ax[1, 0].hist(sel0['a_image'] / sel0['b_image'], 
                              histtype='step', bins=20,
                              color='r')
                ax[1, 1].plot(sel0['a_image'] / sel0['b_image'], 
                              sel0['theta_image'], 'ro')
                ax[0, 0].set_ylabel('N')
                ax[1, 0].set_ylabel('N')        
                ax[1, 0].set_xlabel(r'$\frac{a}{b}$')
                ax[0, 1].set_ylabel(r'$\theta$')
                ax[1, 1].set_ylabel(r'$\theta$')
                ax[1, 1].set_xlabel(r'$\frac{a}{b}$')
                plt.suptitle(r'$\frac{a}{b}$ vs $\theta$ statistics. Cuts applied')
                plt.subplots_adjust(wspace=0.4)
                plt.show()
            # Coordinates of a star I know has horizontal banding
            # x0, y0 = 1095, 1600
            # x0, y0 = 386, 1527 # left amplifier good source
            # x0, y0 = 1253, 3068
            # x0, y0 = 1120, 3554
            if False:
                # Get its params by locating the closest entry in the table
                p = np.power(dfcat['x_image'] - x0, 2)
                p += np.power(dfcat['y_image'] - y0, 2)
                p = np.sqrt(p)
                argmin_p = p.idxmin()
                # Print position arguments
                print(dfcat.iloc[argmin_p][['x_image', 'y_image', 
                                           'xmin_image', 'xmax_image',
                                           'ymin_image', 'ymax_image',
                                           'a_image', 'b_image', 'theta_image',
                                           'number']])
                print('Distance = {0:.3f} pix'.format(p.iloc[argmin_p]))
                # Get needed parameters for displaying the region around the source 
                xcnt, ycnt, xmin, xmax, ymin, ymax, a, b, theta = dfcat.iloc[
                    argmin_p][['x_image', 'y_image',
                    'xmin_image', 'xmax_image', 
                    'ymin_image', 'ymax_image',
                    'a_image', 'b_image', 'theta_image']
                ].values.astype(int)
                stamp = sci[ymin - 1:ymax - 1, xmin - 1:xmax -1]
                xcnt -= xmin
                ycnt -= ymin
                xcnt, ycnt = int(xcnt), int(ycnt)
            # -----------------------------------------------------------------
        t_f1 = time.time()
        print('{0:.2f} sec'.format((t_f1 - t_f0)))
    t1 = time.time()
    print('Elapsed time: {0:.2f} min'.format((t1 - t0) / 60.))
    # Save Dataframe
    df = pd.DataFrame(res, 
                      columns=['expnum', 'ccdnum', 'nite', 'band', 't_eff',
                               'section',
                               'x_image', 'y_image', 'mag_aper_4',
                               'a_image', 'b_image', 'theta_image', 
                               'number',
                               'mean', 'median', 'std', 'var', 'mad', 
                               'ptp', 'npix',
                               'type', 'LR']
    )
    df.to_csv(outname, index=False, header=True)
    print('Saved: {0}'.format(outname))

if __name__ == '__main__':
    print(time.ctime())

    # Here: argparse for a list of immasked and cat files
    hgral = 'Script to evaluate the symmetry in flux for left/right regions'
    hgral += ' of detected objects'
    arg = argparse.ArgumentParser(description=hgral)
    h0 = 'Input a CSV table containing at least columns: t_eff, path, filetype'
    arg.add_argument('--tab', help=h0)
    aux_ccd = [41, 28, 11, 52]
    h1 = 'List of space-separated CCDs for which to do the calculations.'
    h1 += ' Default: {0}'.format(aux_ccd)
    arg.add_argument('--ccd', help=h1, nargs='+', type=int, default=aux_ccd)
    aux_region = [1, 4096, 1024, 1536] # before end was 1948
    h2 = 'Section of each CCD to be used, as space-separated values'
    h2 += ' Format: x0 x1 y0 y1. Where y-axis'
    h2 += ' is the long dimension. Default: {0}'.format(aux_region)
    arg.add_argument('--area', help=h2, type=int, nargs='+',
                     default=aux_region)
    aux_out = '2region_stat_PID{0}.csv'.format(os.getpid())
    h3 = 'Output name for the CSV table of the calculations. Default:'
    h3 += '{0}'.format(aux_out)
    arg.add_argument('--out', help=h3, default=aux_out)
    # Parse
    arg = arg.parse_args()
    # Note the path is not the specific file location but the parent directory
    aux_main(outname=arg.out, 
             d0_range=arg.area[:2], 
             d1_range=arg.area[2:],
             table=arg.tab,
             ccdlist=arg.ccd)

