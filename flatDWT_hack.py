'''flatDWT
Created: September 29, 2016

This script must be able to detect if a flat is acceptable inmediatly after
exposure
Use dflats tagged as bad in FLAT_QA
Use cropped pieces of code from flatStat and decam_test

Oct 4th: DMWY as selected wavelet (symmetric,orthogonal,biorthogonal)

STEPS:
1) Do it on a single ccd with all the possibilities --done
2) do it well on FP with all the possibilities --done
3) Do it for good/bad flats and compare results

HACK version
Date: April 26, 2017
this version was modified to be performed on the simulated dome flat features
'''
import os
import sys
import time
import numpy as np
import scipy.stats
import scipy.signal
import scipy.interpolate
import matplotlib.pyplot as plt
import fitsio
import pywt
import tables

class Toolbox():
    '''methods to be inserted in other
    '''
    @classmethod
    def detect_outlier(cls,imlayer):
        ''' Estimate Iglewicz and Hoaglin criteria for outlier
        and replace by NaN
        http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
        Formula:
        Z=0.6745(x_i - median(x)) / MAD
        if abs(z) > 3.5, x_i is a potential outlier
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
        '''
        #from statsmodels.robust import scale #alternative mad calculation
        from scipy.stats import norm
        '''Percent point function (inverse of cdf -- percentiles) of a normal
        continous random variable
        '''
        cte = norm.ppf(0.75) #aprox 0.6745
        '''Flatten the image, to estimate median
        '''
        flat_im = imlayer.ravel()
        #flat_im=flat_im[flat_im!=-1]#exclude '-1' values
        MAD = np.median( np.abs( flat_im-np.median(flat_im) ) )
        #alternative: scale.mad(flat_im, c=1, axis=0, center=np.median)
        Zscore = cte*(flat_im-np.median(flat_im))/MAD
        Zscore = np.abs(Zscore)
        ''' search for outliers and if present replace by -1
        '''
        imlayer[np.where(np.abs(cte*(imlayer-np.median(flat_im))/MAD)>3.5)]=-1.
        return imlayer
        '''
        if len(Zscore[Zscore>3.5])>0:
            for k in range(0,imlayer.shape[0]):
                for m in range(0,imlayer.shape[1]):
                    if np.abs( cte*(imlayer[k,m]-np.median(flat_im))/MAD )>3.5:
                        imlayer[k,m] = np.nan
            return imlayer
        else:
            return imlayer
        '''

    @classmethod
    def quick_stat(cls,arr_like):
        MAD = np.median( np.abs(arr_like-np.median(arr_like)) )
        print '__________'
        print '* Min | Max | Mean = {0} | {1} | {2}'.format(
            np.min(arr_like),np.max(arr_like),np.mean(arr_like))
        print '* Median | Std | MAD = {0} | {1} | {2}'.format(
            np.median(arr_like),np.std(arr_like),MAD)
        print '* .25 | .5 | .75 = {0} | {1} | {2}'.format(
            np.percentile(arr_like,.25),np.percentile(arr_like,.5),
            np.percentile(arr_like,.75))
        return False

    @classmethod
    def dwt_library(cls,img_arr):
        '''Display all the availbale DWT (single level) for the input array
        '''
        t1 = time.time()
        count = 0
        for fam in pywt.families():
            for mothwv in pywt.wavelist(fam):
                for mod in pywt.Modes.modes:
                    print '\tWavelet: {0} / Mode: {1}'.format(mothwv,mod)
                    (c_A,(c_H,c_V,c_D)) = pywt.dwt2(img_arr,
                                                    pywt.Wavelet(mothwv),
                                                    mod)
                    count += 1
                    fig=plt.figure(figsize=(6,11))
                    ax1=fig.add_subplot(221)
                    ax2=fig.add_subplot(222)
                    ax3=fig.add_subplot(223)
                    ax4=fig.add_subplot(224)
                    ax1.imshow(c_A,cmap='seismic')#'flag'
                    ax2.imshow(c_V,cmap='seismic')
                    ax3.imshow(c_H,cmap='seismic')
                    ax4.imshow(c_D,cmap='seismic')
                    plt.title('Wavelet: {0} / Mode: {1}'.format(mothwv,mod))
                    plt.subplots_adjust(left=.06, bottom=0.05,
                                        right=0.99, top=0.99,
                                        wspace=0., hspace=0.)
                    plt.show()
        t2 = time.time()
        print ('\nTotal time in all {1} modes+mother wavelets: {0:.2f}\''
               .format((t2-t1)/60.,count))

    @classmethod
    def range_str(cls,head_rng):
        head_rng = head_rng.strip('[').strip(']').replace(':',',').split(',')
        return map(lambda x: int(x)-1, head_rng)

    @classmethod
    def plot_flat(cls,img):
        '''out of the FP points are masked with value -1
        '''
        img = np.ma.masked_where((img==-1),img)
        percn = [np.percentile(img,k) for k in range(1,101,1)]
        ind_aux = len(percn) - percn[::-1].index(-1)
        minVal = percn[ind_aux]
        #
        fig = plt.figure(figsize=(12,6))

        ax1 = plt.subplot2grid((2,3),(0,0),colspan=1,rowspan=1)
        ax2 = plt.subplot2grid((2,3),(0,1),colspan=2,rowspan=2)
        ax3 = plt.subplot2grid((2,3),(1,0),colspan=1,rowspan=1)
        #ax1 = fig.add_subplot(221,adjustable='box',aspect=1.)
        #ax2 = fig.add_subplot(222,adjustable='box',aspect=1.)
        #ax3 = fig.add_subplot(223)
        #contour filled
        from matplotlib import colors,ticker,cm
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        Y = np.linspace(0,img.shape[0]-1,img.shape[0])
        X = np.linspace(0,img.shape[1]-1,img.shape[1])
        X,Y = np.meshgrid(X,Y)
        cs = ax1.contourf(X,Y,img,locator=ticker.LogLocator(),cmap=cm.PuBu_r,
                        nchunk=0,lines='solid',origin='lower')
        divider1 = make_axes_locatable(ax1)
        cax1 = divider1.append_axes('left',size='10%',pad=0.9)
        cbar1 = plt.colorbar(cs,cax=cax1)
        #imshow
        img_aux = img #np.ma.masked_where((img<1),img)
        im = ax2.imshow(img_aux,cmap=cm.PuBu_r,vmin=minVal,origin='lower')
        ax2.contour(X,Y,img,locator=ticker.LogLocator(),colors='red',
                    nchunk=0,linestyles='solid',linewidths=1,origin='lower')
        divider2 = make_axes_locatable(ax2)
        cax2 = divider2.append_axes('right',size='10%',pad=0.2)
        cbar2 = plt.colorbar(im,cax=cax2)
        #histogram and univariate
        imgH = np.sort(img.ravel())
        p,x = np.histogram(imgH[np.where(imgH>minVal)],
                        bins=100)
        #convert bin edges to centers
        x = x[:-1] + (x[1]+x[2])/2.
        f = scipy.interpolate.UnivariateSpline(x,p,s=len(imgH)/20)
        ax3.hist(imgH[np.where(imgH>minVal)],bins=100,
                histtype='stepfilled',color='#43C6DB')
        #
        ax2.invert_yaxis
        ax3.plot(x,f(x),'k-')
        ax3.set_ylim([0,np.max(f(x))])
        plt.tight_layout()
        #plt.subplots_adjust(left=.025,bottom=0.05,right=0.96,top=0.97)
        plt.show()


class FPBinned():
    def __init__(self,folder,fits):
        '''Simple method to open focal plane binned images
        When a position not belongs to focal plane, the value is -1
        Before return it, add 1 to set outer region to zero value
        '''
        fname = os.path.join(folder,fits)
        tmp = np.load(fname)
        #outliers must be present, because there are not removed until
        #later steps, tmp = Toolbox.detect_outlier(tmp)
        tmp += 1.
        self.fpBinned = tmp


class DWT():
    @classmethod
    def single_level(cls,img_arr,wvfunction='dmey',wvmode='zero'):
        '''DISCRETE wavelet transform
        Wavelet families available: 'haar', 'db', 'sym', 'coif', 'bior',
        'rbio','dmey'
        http://www.pybytes.com/pywavelets/regression/wavelet.html
        - When flat shows issues, it presents discontinuities in flux
        - To perform the wavelet, border effects must be considered
        - Bumps at the edge are by now flagged but if something
          could be done, it would represent a huge improvement, specially
          for th SN team
        FOR THE ENTIRE FOCAL PLANE
        --------------------------
        Must fill the interCCD space with zeroes (?) or interpolation.
        Scales must define the refinement scales I will look. The more scales,
        the slower calculation and better the resolution.
        - MODES: different ways to deal with border effects. Zeros is preferred
        - DWT2/WAVEDEC2 output is a tuple (cA, (cH, cV, cD)) where (cH, cV, cD)
          repeats Nwv times
        Coeffs:
        c_A : approximation (mean of coeffs) coefs
        c_H,c_V,c_D : horizontal detail,vertical,and diagonal coeffs
        '''
        (c_A,(c_H,c_V,c_D)) = pywt.dwt2(img_arr,pywt.Wavelet(wvfunction),
                                        wvmode)
        '''reduced set of parameters through pywt.threshold()
        with args: soft, hard, greater, less
        To define a set of points (with same dimension as the image), the
        detailed coeffs must be employed
        '''
        return c_A,c_H,c_V,c_D

    @classmethod
    def dec_mlevel(cls,img_arr,wvfunction='dmey',wvmode='zero',Nlev=8):
        '''Wavelet Decomposition in multiple levels (decimated), differs to
        DWT which is the wavelet transform of one level only
        Inputs:
        - img_arr: 2D array containing image data
        - wvfunction: mother wavelet to be used
        - wvmode: method to deal with borders
        - Nlev: number of level for decomposition (max is 8)
        Outputs:
        - WAVEDEC2 output is a tuple (cA, (cH, cV, cD)) where (cH, cV, cD)
          repeats Nwv times
        '''
        c_ml = pywt.wavedec2(img_arr,pywt.Wavelet(wvfunction),
                             wvmode,level=Nlev)
        aux_shape = []
        for i in range(len(c_ml)):
            if (i == 0): aux_shape.append(c_ml[i].shape)
            else: aux_shape.append(c_ml[i][0].shape)
        cls.cmlshape = aux_shape
        return c_ml

    @classmethod
    def undec_mlevel(cls,img_arr,wvfunction='dmey',Nlev=8):
        '''Stationary wavelet decomposition in 2D, also
        '''
        c_ml = pywt.swt2(img_arr,pywt.Wavelet(wvfunction),Nlev,
                        start_level=0,axes=(0,1))
        print len(c_ml)
        return c_ml


class Coeff(DWT):
    '''method for save results of DWT on a compressed pytables
    '''
    @classmethod
    def set_table(cls,tablename,Nlev):
        '''Method for initialize the file to be filled with the results from
        the DWT decomposition. It works either for 2 or 8 levels of
        decomposition.
        Descriptions for pytables taken from "Numerical Python: A practical
        techniques approach for Industry"
        '''
        #create a new pytable HDF5 file handle. This does not represents the
        #root group. To access the root node must use cls.h5file.root
        cls.h5file = tables.open_file(
            tablename,mode='w',
            title='HDF5 file containing data of dome flat wavelet-behavior.',
            driver='H5FD_CORE')
        #create groups of the file handle object. Args are: path to the parent
        #group (/), the group name, and optionally a title for the group. The
        #last is a descriptive attribute can be set to the group.
        group = cls.h5file.create_group(
            '/','dwt',
            title='Created using PyWavelets {0}'.format(pywt.__version__))
        #the file handle is defined, also the group inside it. Under the group
        #we will save the DWT tables.
        #to create a table with mixed type structure we create a class that
        #inherits from tables.IsDescription class
        if Nlev == 2:
            class Levels(tables.IsDescription):
                c_A = tables.Float32Col(shape=DWT.cmlshape[0])
                c1 = tables.Float32Col(shape=DWT.cmlshape[1])
                c2 = tables.Float32Col(shape=DWT.cmlshape[2])
        elif Nlev == 8:
            class Levels(tables.IsDescription):
                c_A = tables.Float32Col(shape=DWT.cmlshape[0])
                c1 = tables.Float32Col(shape=DWT.cmlshape[1])
                c2 = tables.Float32Col(shape=DWT.cmlshape[2])
                c3 = tables.Float32Col(shape=DWT.cmlshape[3])
                c4 = tables.Float32Col(shape=DWT.cmlshape[4])
                c5 = tables.Float32Col(shape=DWT.cmlshape[5])
                c6 = tables.Float32Col(shape=DWT.cmlshape[6])
                c7 = tables.Float32Col(shape=DWT.cmlshape[7])
                c8 = tables.Float32Col(shape=DWT.cmlshape[8])
        else:
            raise ValueError('Method supports 2 or 8 decomposition levels.')
        #with the table structure already defined, the table with DWT results
        #can be fully created. Args are: a group object or the path to the root
        #node, the table name, the table structure specification, and
        #optionally the table title. The last is stored as attribute.
        cls.cml_table = cls.h5file.create_table(
            group,'dmeyN2',Levels,title='Mother wavelet:dmey. Levels:2.')

    @classmethod
    def fill_table(cls,coeff_tuple,Nlev,dict_database):
        '''Method to fill the HDF5 file with the DWT results. It works for 2
        and 8 levels of decomposition.
        '''
        #using .attrs or ._v_attrs we can access different levels in the
        #HDF5 structure
        cls.cml_table.attrs.DB_INFO = dict_database
        cml_row = Coeff.cml_table.row
        if Nlev == 2:
            for m in xrange(3):
                cml_row['c_A'] = coeff_tuple[0]
                cml_row['c1'] = coeff_tuple[1][m]
                cml_row['c2'] = coeff_tuple[2][m]
                cml_row.append()
                Coeff.cml_table.flush()
        elif Nlev == 8:
            for m in xrange(3):
                cml_row['c_A'] = coeff_tuple[0]
                cml_row['c1'] = coeff_tuple[1][m]
                cml_row['c2'] = coeff_tuple[2][m]
                cml_row['c3'] = coeff_tuple[3][m]
                cml_row['c4'] = coeff_tuple[4][m]
                cml_row['c5'] = coeff_tuple[5][m]
                cml_row['c6'] = coeff_tuple[6][m]
                cml_row['c7'] = coeff_tuple[7][m]
                cml_row['c8'] = coeff_tuple[8][m]
                cml_row.append()
                Coeff.cml_table.flush()
        else:
            raise ValueError('Method supports 2 or 8 decomposition levels.')

    @classmethod
    def close_table(cls):
        Coeff.cml_table.flush()
        Coeff.h5file.close()


if __name__=='__main__':
    #for binning of 64x64 the maximum level in wavedec2 is 2
    #for binning of 32x32 the maximum level in wavedec2 is 3
    #for binning of 16x16 the maximum level in wavedec2 is 4
    #for binning of 8x8 the maximum level in wavedec2 is 5
    #for binning of 4x4 the maximum level in wavedec2 is 6
    #testing...
    path = os.path.join(os.path.expanduser('~'),
                    'Code/des_calibrations/dwt_test_files')
    fn = lambda s: os.path.join(path,s)
    bin_fp = np.load(fn('feat2_b44.npy'))
    print bin_fp.shape
    c1_ml = DWT.dec_mlevel(bin_fp,Nlev=8)
    c2_ml = DWT.undec_mlevel(bin_fp)

    '''Coeff.set_table(fnout,decLev)
    #fill table
    Coeff.fill_table(c_ml,decLev,binfo)
    #close table
    Coeff.close_table()
    count += 1
    t2 = time.time()
    if not count%10:
        print '\n{2}--{1},read+multilevel+store: {0:.2f}\'\''.format(
            t2-t1,dflat_tab['filename'][k],count)
    '''
