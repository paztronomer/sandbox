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
import logging
import socket
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
                logging.info('Creating directory {0}'.format(folder))
            except:
                logging.info("Issue creating the folder {0}".format(folder))
        else:
            logging.info('Checking: directory {0} exists'.format(folder))
        return True

    @classmethod
    def split_path(cls,path):
        '''Method to return the relative root folder (one level upper),
        given a path.
        '''
        #relat_root = os.path.abspath(os.path.join(path,os.pardir))
        relroot,filename = os.path.split(path)
        return relroot,filename


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
    def old_set_table(cls,dwt_res,table_name=None,guess_rows=8,
                    label='',title=''):
        '''Method for initialize the file to be filled with the results from
        the DWT decomposition.
        Descriptions for pytables taken from "Numerical Python: A practical
        techniques approach for Industry"
        '''
        if table_name is None:
            table_name = 'pid{0}.table'.format(os.getpid())
        #create a new pytable HDF5 file handle. This does not represents the
        #root group. To access the root node must use cls.h5file.root
        cls.h5file = tables.open_file(
            table_name,
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
        #If use a ndarray, the created table only assumes the data type, not
        #the shape
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
            title=title,
            expectedrows=guess_rows)
        #cls.h5file.flush()
        #https://groups.google.com/forum/#!topic/pytables-users/EqxD5zHroc8

    @classmethod
    def set_table(cls,dwt_res,table_name=None,guess_rows=8,label='',
                    title=''):
        '''Method for initialize the file to be filled with the results from
        the DWT decomposition.
        Descriptions for pytables taken from "Numerical Python: A practical
        techniques approach for Industry"
        '''
        if table_name is None:
            table_name = 'pid{0}.table'.format(os.getpid())
        #create a new pytable HDF5 file handle. This does not represents the
        #root group. To access the root node must use cls.h5file.root
        cls.h5file = tables.open_file(
            table_name,
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

        #Will use a dictionary, to setup columns in native Col()
        #Mind the structure: [(cA,(cH,cV,cD)),(cA,(cH,cV,cD)),...]
        diccio = dict()
        for i in xrange(len(dwt_res)):
            coeff = 'coeff{0}'.format(i+1)
            diccio[coeff] = tables.FloatCol(shape=dwt_res[i][0].shape)

        #with the table structure already defined, the table with DWT results
        #can be fully created. Args are: a group object or the path to the root
        #node, the table name, the table structure specification, and
        #optionally the table title. The last is stored as attribute.
        cls.cml_table = cls.h5file.create_table(
            group,
            label,
            description=diccio,
            title=title,
            expectedrows=guess_rows)
        #cls.h5file.flush()
        #https://groups.google.com/forum/#!topic/pytables-users/EqxD5zHroc8

    @classmethod
    def fill_table(cls,dwt_res,info_table={}):
        '''Method to fill the HDF5 file with the DWT results
        '''
        #using .attrs or ._v_attrs we can access different levels in the
        #HDF5 structure
        cls.cml_table.attrs.DB_INFO = info_table
        cml_row = Coeff.cml_table.row
        #fill row by row, using the first level as template. Beforehand we
        #know every level of decmposition will has the following setup:
        #(cA,(cH,cV,cD))
        for lev in xrange(len(dwt_res)):
            decomp = dwt_res[lev]
            cml_row['coeff{0}'.format(lev+1)] = decomp[0]
            cml_row.append()
            cml_row['coeff{0}'.format(lev+1)] = decomp[1][0]
            cml_row.append()
            cml_row['coeff{0}'.format(lev+1)] = decomp[1][1]
            cml_row.append()
            cml_row['coeff{0}'.format(lev+1)] = decomp[1][2]
            cml_row.append()
            Coeff.cml_table.flush()

    @classmethod
    def close_table(cls):
        Coeff.cml_table.flush()
        Coeff.h5file.close()

class Caller(Coeff):
    """Class inherites from Coeff, which inherites from DWT
    """
    def __init__(self,fname,wvmother):
        t0 =  time.time()
        root,npyfile = Toolbox.split_path(fname)
        info = '\nWorking on {0}. {1}'.format(npyfile,time.ctime())
        info += '\n\t(*)Wavelet: {0}'.format(wvmother)
        print info
        logging.info(info)
        #load the binned focal plane
        bin_fp = np.load(fname)
        #for each mother wavelet, save in a different folder
        out_folder = os.path.join(root,wvmother)
        #check if the folder exists, if not, crete it
        Toolbox.check_folder(out_folder)
        #setup the output name
        out_name = '{0}_{1}.h5'.format(wvmother,npyfile[:npyfile.find('.npy')])
        out_name = os.path.join(out_folder,out_name)
        #perform the DWT
        res = DWT.undec_mlevel(bin_fp,wvfunction=wvmother)
        #setup the table to harbor results
        Coeff.set_table(res,
                    table_name=out_name,
                    label='fpazch_at_illinois_dot_edu',
                    title='DWT type: {0}'.format(wvmother))
        #fill the table
        d = dict()
        d['module'] = 'PyWavelets v{0}'.format(pywt.__version__)
        d['wavelet'] = wvmother
        d['time'] = time.ctime()
        Coeff.fill_table(res,info_table=d)
        #close the file
        Coeff.close_table()
        t1 = time.time()
        info_end = 'Ended with {0}. Elapsed: {1:.3f} sec'.format(npyfile,t1-t0)
        print info_end
        logging.info(info_end)


if __name__=='__main__':
    #directory setup
    #campus cluster precal nodes or macbook
    if socket.gethostname() in ['ccc0027','ccc0028','ccc0029']:
        npy_folder = 'scratch/dwt_files'
    elif  socket.gethostname() in ['albatross.local']:
        npy_folder = 'Code/des_calibrations/dwt_files'
    else:
        logging.warning('Need to define the location of Binned FP')

    path = os.path.join(os.path.expanduser('~'),npy_folder)
    fn = lambda s: os.path.join(path,s)
    #wavelet library setup
    #wvmother = ['dmey','morl','cmor','shan','gaus','mexh','haar']
    wvmother = ['dmey','haar']

    for root,dirs,files in os.walk(path):
        FPname = [fn(binned) for binned in files if ('.npy' in binned)]
        '''below can be replaced by a map()'''
        for fn in FPname:
            for wv_n in wvmother:
                if len(pywt.wavelist(wv_n)) == 1:
                    doit = Caller(fn,wv_n)
                elif len(pywt.wavelist(wv_n)) > 1:
                    for wv_sub in pywt.wavelist(wv_n):
                        doit = Caller(fn,wv_sub)


    exit(0)

    '''For testing on the DWT perform'''
    #for binning of 64x64 the maximum level in wavedec2 is 2
    #for binning of 32x32 the maximum level in wavedec2 is 3
    #for binning of 16x16 the maximum level in wavedec2 is 4
    #for binning of 8x8 the maximum level in wavedec2 is 5
    #for binning of 4x4 the maximum level in wavedec2 is 6

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

    Coeff.set_table(coeff,table_name='test.table',label='b3232',
                title='DWT type: {0}'.format(wvmother))
    '''fill the table'''
    d = dict()
    d['module'] = 'PyWavelets v{0}'.format(pywt.__version__)
    d['wavelet'] = wvmother
    d['time'] = time.ctime()
    Coeff.fill_table(coeff,info_table=d)
    Coeff.close_table()
