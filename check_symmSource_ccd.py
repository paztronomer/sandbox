""" Script to evaluate the circular symmetry of sources 
"""

import os
import glob
import time
import argparse
import logging
import numpy as np
from scipy import interpolate
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
              angle_lr=[np.pi, 0], 
              alpha=10, 
              sampling_factor=3,
              img=None):
    """ Method to get data from 2 circular regions at the left and right
    or a centroid. Returns 2 masked arrays, one for the left and other for
    the right side of the centroid.
    The values of the fine grid images have its values normalized to sum the
    initial value of its parent pixel.
    This method works well if:
    * borders of circular regions are not +/- pi/2
    * both regions can be defined with only 2 straight lines
    Inputs
    - cnt_coo: centroid, where coordinates are [x_image, y_image]
    - radius: adjusted radio of the object
    - angle_lr: left, right reference angles to create the circular regions 
    - alpha: angle in degrees to set a range above and below the left/right
    reference angles
    """
    # Get the circle sectors, for left and right of the centroid
    # Radius for statistics
    radius = 2 * a
    # Circular section 
    # Upper and lower angles
    theta_l_up = angle_lr[0] - np.deg2rad(alpha)
    theta_l_dw = angle_lr[0] + np.deg2rad(alpha)
    theta_r_up = angle_lr[1] + np.deg2rad(alpha)
    theta_r_dw = angle_lr[1] - np.deg2rad(alpha)
    # (x, y) coordinates of limit positions, from polar coordinates
    xl = radius * np.cos(theta_l_up)
    yl_up, yl_dw = radius * np.sin(theta_l_up), radius * np.sin(theta_l_dw)
    xr = radius * np.cos(theta_r_up)
    yr_up, yr_dw = radius * np.sin(theta_r_up), radius * np.sin(theta_r_dw)
    # Round values for coordinates limits, as they represent pixel positions
    # xl, yl_up, yl_dw = np.round(xl), np.round(yl_up), np.round(yl_dw)
    # xr, yr_up, yr_dw = np.round(xr), np.round(yr_up), np.round(yr_dw)
    # xl, yl_up, yl_dw = int(xl), int(yl_up), int(yl_dw)
    # xr, yr_up, yr_dw = int(xr), int(yr_up), int(yr_dw)
    #
    # Discretize just at the very end
    #
    # Translate the coordinates, by the offset of the centroid
    xcnt, ycnt = cnt_coo
    xl += xcnt
    yl_up += ycnt
    yl_dw += ycnt 
    xr += xcnt
    yr_up += ycnt
    yr_dw += ycnt
    # The precission value 0.001 was defined arbitrary
    if (((yl_dw - yr_dw) < 0.001) and ((yl_up - yr_up) < 0.001)):
        pass
    else:
        logging.error('Circular sections are not symmetric')
        exit(1)
    # Get the pixels inside the circle sections: 2 lines needs to be adjusted
    # Outputs from linregress are: m, b, r_value, p_value, std_err
    f = interpolate.interp1d([xl, xcnt, xr], 
                             [yl_up, ycnt, yr_dw], 
                             kind='linear',
                             bounds_error=False,
                             fill_value='extrapolate')
    g = interpolate.interp1d([xl, xcnt, xr], 
                             [yl_dw, ycnt, yr_up], 
                             kind='linear',
                             bounds_error=False,
                             fill_value='extrapolate')
    # Evaluate each one of the coordinates to see if they belong to circle 
    # left/right sections
    #
    # Resampling the grid to make it finer
    #
    # Increase the number of points of the grid to each pixel be 3x3 or
    # the value defined by the sampling_factor
    xmin, xmax = int(np.floor(xl)), int(np.ceil(xr))
    ymin, ymax = int(np.floor(yl_dw)), int(np.ceil(yl_up))
    Nx = (sampling_factor) * (xmax - xmin) + 1
    Ny = (sampling_factor) * (ymax - ymin) + 1
    # Define the supersampled grid
    yy, xx = np.meshgrid(np.linspace(ymin, ymax, Ny), 
                         np.linspace(xmin, xmax, Nx),
                         sparse=False, indexing='ij')
    #
    # Interpolate grid values
    #
    # 1) Define initial grid from the image, and transform it to set of 
    # coordinates
    yy_ini, xx_ini = np.mgrid[ymin:ymax, xmin:xmax]
    points = np.vstack([xx_ini.ravel(), yy_ini.ravel()])
    values = img[ymin:ymax, xmin:xmax].flatten()
    # 2) Using the positions and values from the original image, flesh out the
    # finer grid. Nearest value will be assumed
    interp_img = interpolate.griddata(points.T, values, (xx, yy), 
                                      method='nearest')
    # Normalize byt the number of subpixels each pixel is divided
    interp_img /= np.power(sampling_factor, 2)
    # Plotting the resampled grid for evaluation
    if not True:
        im_norm = ImageNormalize(img, 
                                 interval=ZScaleInterval(),
                                 stretch=SqrtStretch())
        im_norm2 = ImageNormalize(interp_img, 
                                  interval=ZScaleInterval(),
                                  stretch=SqrtStretch())
        fig, ax = plt.subplots(1, 3, figsize=(9, 3))
        ax[0].imshow(img, norm=im_norm, origin='lower')
        ax[1].imshow(img, norm=im_norm, origin='lower', cmap='viridis')
        ax[1].scatter(xx, yy, marker='.', s=5, color='white', alpha=0.3)
        ax[2].imshow(interp_img, norm=im_norm2, origin='lower')
        #
        ax[0].set_xlim([xmin, xmax])
        ax[0].set_ylim([ymin - 10, ymax + 10])
        ax[1].set_xlim([xmin, xmax])
        ax[1].set_ylim([ymin - 10, ymax + 10])
        plt.show()
    #
    # Select the points from the fine grid belonging to the region of interest 
    # Create a mask
    #
    c_left = np.logical_and(yy <= f(xx), yy >= g(xx))
    c_right = np.logical_and(yy <= g(xx), yy >= f(xx))
    c_both = np.logical_or(c_left, c_right)
    # Mask image
    ma_img = np.ma.masked_where(~c_both, interp_img)
    ma_left = np.ma.masked_where(~c_left, interp_img)
    ma_right = np.ma.masked_where(~c_right, interp_img)
    # Return only left/right
    return ma_left, ma_right

def masked_stat(ma_x):
    """ Receives a masked array and performs some basic statistics
    """
    ma_x = ma_x.compressed()
    mad = lambda x: np.median(np.abs(x - np.median(x)))
    out = [mad(ma_x), np.mean(ma_x), np.std(ma_x), ma_x.size]
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
    # Flatten the arrays, because glob.glob introducesan additional level of
    # nesting
    fnm1 = np.array(fnm1).flatten()
    fnm2 = np.array(fnm2).flatten()

    # Results list to construct the dataframe
    res = []

    for i in range(len(fnm1)):
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
        # Alse save T_EFF 
        t_eff = df.loc[df['expnum'] == expnum, 't_eff'].values[0]

        for obj in range(cnt_x.size): 
            # Using the above parameters define 2 circular regions 
            try:
                r_left, r_right = region_lr([cnt_x[obj], cnt_y[obj]], a[obj], 
                                            angle_lr=[np.pi, 0], 
                                            alpha=22.5, 
                                            img=sci)
            except:
                print(sel.iloc[obj]['flags'])
                exit()
            # Get the statistics for each region
            st_l = masked_stat(r_left)
            st_r = masked_stat(r_right)
            # Get fitted parameters for each object
            fpar = [cnt_x[obj], cnt_y[obj], a[obj], 
                    sel['b_image'].values[obj],
                    sel['theta_image'].values[obj], 
                    sel['number'].values[obj]]
            #
            aux_l = [expnum, ccdnum, nite, band, t_eff, region] + fpar + st_l
            aux_l += ['L']
            aux_r = [expnum, ccdnum, nite, band, t_eff, region] + fpar + st_r
            aux_r += ['R']
            #
            res.append(aux_l)
            res.append(aux_r)

            # Seems that a higher MAD and lower mean is a good test for a healthy
            # stellar profile
            # Run on a set of bad and good exposures.
            
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
    t1 = time.time()
    print('Elapsed time: {0:.2f} min'.format((t1 - t0) / 60.))
    # Save Dataframe
    df = pd.DataFrame(res, 
                      columns=['expnum', 'ccdnum', 'nite', 'band', 't_eff',
                               'x_image', 'y_image', 
                               'a_image', 'b_image', 'theta_image', 
                               'number',
                               'region', 'mad', 'mean', 'std', 'npix',
                               'LR']
    )
    # Need to include centroid, a, b, theta, object ID
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
    aux_ccd = [41, 28, 4, 55]
    h1 = 'List of space-separated CCDs for which to do the calculations.'
    h1 += ' Default: {0}'.format(aux_ccd)
    arg.add_argument('--ccd', help=h1, nargs='+', type=int, default=aux_ccd)
    aux_region = [1, 4096, 1024, 1948]
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

