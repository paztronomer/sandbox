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
import pickle
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
    def range_str(cls,head_rng):
        head_rng = head_rng.strip('[').strip(']').replace(':',',').split(',')
        return map(lambda x: int(x)-1, head_rng)

    @classmethod
    def check_folder(cls,folder):
        '''Method to check for the folder existence, if not present, tries to
        crete it
        '''
        if not os.path.exists(folder):
            try:
                os.makedirs(folder)
            except:
                logging("Issue cueatring the folder {0}".format(folder))
        return True

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
        while True:
            try:
                c_ml = pywt.wavedec2(img_arr,pywt.Wavelet(wvfunction),
                                    wvmode,level=Nlev)
                break
            except:
                Nlev -= 1
        return c_ml

    @classmethod
    def undec_mlevel(cls,img_arr,wvfunction='dmey',lev_end=8):
        '''Stationary wavelet decomposition in 2D. Undecimated.
        Inputs:
        - img_arr: 2D array
        - wvfunction: mother wavelet, string
        - lev_end: final enad of decomposition. Maximum is 8
        - lev_ini: initial level of decomposition. Minumum is 0
        Outputs:
        -
        '''
        while True:
            try:
                c_ml = pywt.swt2(img_arr,pywt.Wavelet(wvfunction),
                                lev_end,axes=(0,1))
                break
            except:
                lev_end -= 1
        return c_ml


class Coeff(DWT):
    '''method for save results of DWT on a compressed pytables
    '''
    @classmethod
    def set_table(cls,dwt_res,tablename,guess_rows=4,label='',title=''):
        '''Method for initialize the file to be filled with the results from
        the DWT decomposition.
        Descriptions for pytables taken from "Numerical Python: A practical
        techniques approach for Industry"
        '''
        #create a new pytable HDF5 file handle. This does not represents the
        #root group. To access the root node must use cls.h5file.root
        cls.h5file = tables.open_file(
            tablename,
            mode='w',
            title='HDF5 table for DWT',
            driver='H5FD_CORE')
        #create groups of the file handle object. Args are: path to the parent
        #group (/), the group name, and optionally a title for the group. The
        #last is a descriptive attribute can be set to the group.
        group = cls.h5file.create_group(
            '/',
            'dflat',
            title='Dome Flat Analysis')
        #the file handle is defined, also the group inside it. Under the group
        #we will save the DWT tables.

        #to create a table with mixed type structure we create a class that
        #inherits from tables.IsDescription class, or we can create from a
        #dictionary, ndarray, among others.
        #See http://www.pytables.org/usersguide/libref/
        #file_class.html#tables.File.create_table
        #To initialize the table, will use a structurred array
        #Mind the structure: [(cA,(cH,cV,cD)),(cA,(cH,cV,cD)),...]
        dt = []
        for idx,i in enumerate(dwt_res):
            dt.append(('coeff{0}'.format(idx+1),'f8'))
            if idx == 0:
                aux = i[0]
            else:
                aux = np.dstack((aux,i[0]))
        aux = np.array(aux,dtype=dt)
        #with the table structure already defined, the table with DWT results
        #can be fully created. Args are: a group object or the path to the root
        #node, the table name, the table structure specification, and
        #optionally the table title. The last is stored as attribute.
        cls.cml_table = cls.h5file.create_table(
            group,
            label,
            description=aux,
            title=title)
        #cls.h5file.flush()
        #https://groups.google.com/forum/#!topic/pytables-users/EqxD5zHroc8

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
    path = os.path.join(os.path.expanduser('~'),
                    'Code/des_calibrations/dwt_test_files')
    fn = lambda s: os.path.join(path,s)
    bin_fp = np.load(fn('feat2_b3232.npy'))

    '''For testing on the DWT perform'''
    #for binning of 64x64 the maximum level in wavedec2 is 2
    #for binning of 32x32 the maximum level in wavedec2 is 3
    #for binning of 16x16 the maximum level in wavedec2 is 4
    #for binning of 8x8 the maximum level in wavedec2 is 5
    #for binning of 4x4 the maximum level in wavedec2 is 6
    #testing...
    wvmother = 'dmey'
    """
    t1 =  time.time()
    c1_ml = DWT.dec_mlevel(bin_fp)
    print time.time()-t1
    t2 = time.time()
    c2_ml = DWT.undec_mlevel(bin_fp,wvfunction=wvmother)
    print time.time()-t2
    """
    #pickle.dump(c2_ml,open('coeff.pickle','w+'))
    coeff = pickle.load(open('coeff.pickle','r'))

    '''To test on the table creation'''
    Coeff.set_table(coeff,'test.table',label='b3232',
                title='DWT type: {0}'.format(wvmother))

    exit()

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
