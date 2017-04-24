'''
Script to create mimic dome flats, by seeting up a focal plane array with
some feature and then applying the mask of the CCD sections, which effectively
catches the light
Different to my usual codes, I'll use self instead of cls
'''
import os
import sys
import logging
import time
import numpy as np
import scipy
import scipy.interpolate
#from scipy import ndimage
import tables
import fitsio
import pickle
#setup for display visualization
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

'''
STEPS
1) read and pickle header from an exposure, to get the CCDs positions --DONE
2) mask the borders of each CCD --DONE
3) create features
4) resize Fp and mask to different scales
Pending:
1) crop the borders, np.all(a == a[0,:], axis = 0)
'''


class CCD_position():
    def __init__(self,path='/work/devel/fpazch/calib_space',
                fname='DECam_00615355.fits.fz',
                rerun=False):
        '''Reads the header of a raw exposure and gets the information for
        CCD positioning
        Extensions runs from 0 to 70, where 0 is the non-data
        default extension. The division is:
        - 1 to 62: science / hduname '{N/S}{1-31}'
        - 63 to 70: guider/focus / hduname 'F{N/S}{1-4}'

        To read data: hdu[ext][:,:] or hdu[ext].read()
        To read header keys: hdu[ext].read_header()['key']
        '''
        o1,o2 = 'ccd_dim.pickle','fp_dim.pickle'
        if rerun:
            d = os.path.join(path,fname)
            M_hdu = fitsio.FITS(d)
            dim = np.empty((0,5))
            for M in M_hdu[1:]:
                x = M.read_header()
                if (x['CCDNUM'] <= 62) and (x['CCDNUM'] > 0):
                    x0 = CCD_position().dim_str(x['DETSIZE'])
                    x1 = np.array([x['CCDNUM']])
                    x2 = CCD_position().dim_str(x['DETSEC'])
                    line = np.concatenate([x1,x2])[np.newaxis,:]
                    dim = np.vstack([dim,line])
                else:
                    pass
                    #63 to 70 are guider/focus
            fp_dim = np.r_[x0[1],x0[3]]
            dim = dim.astype(int)
            pickle.dump(dim,open(o1,'w+'))
            pickle.dump(fp_dim,open(o2,'w+'))
            self.ccd_dim = dim
            self.tot_dim = fp_dim
        else:
            self.ccd_dim = pickle.load(open(o1,'r'))
            self.tot_dim = pickle.load(open(o2,'r'))

    def dim_str(self,str_range):
        res = str_range.strip('[').strip(']').replace(':',',').split(',')
        return np.array(map(int,res))


class Mimic(): #(CCD_position):
    def crop_borders(self,arr):
        '''Given a 2D array, returns the lowest/highest coordinates for which
        all elements are zeros, so we can reduce the dimension by throwing
        away the unused borders. This method must check no non-zero pixels
        apre present above the threshold
        '''
        #r[~np.all(r == 0, axis=1)]
        #r[~np.all(r == 0, axis=1)]
        return False

    def feature_fp(self):
        '''Method to create non uniformity shapes in the FP, to test the
        detection efficiency of the DWT algorithm at different scales
        '''
        size = 0
        #np.random.multivariate_normal(mean, cov, 10000)
        return False

    def join_binned(self,data=None,mask=None,binsize=16,header_again=False):
        """Method to reduce the dimensions of the simulated focal plane and
        of the mask. To not mess with combining data from mask and image,
        both of them will be processed separately and then joined
        Inputs:
        - data
        - mask
        - binsize: if is an integer, use square bin, if a tuple, use its
        dimension
        """
        if isinstance(binsize,(int,long)):
            binsize = (binsize,binsize)
        elif (isinstance(binsize,tuple)) and (len(binsize)==2):
            binsize = binsize
        else:
            logging.warning("Binning: dim must be either 1 or 2 integers")
            exit(1)

        #scale mask
        #as mask is a binary-valued object, nearest neighbor interpolation is
        #the selected algorithm
        '''put this in the caller'''
        mask = Mimic().ccd_mask()
        x0 = np.floor(mask.shape[0]/np.float(binsize[0])).astype(int)
        x1 = np.floor(mask.shape[1]/np.float(binsize[1])).astype(int)
        s_mask = Mimic().scale_mask(mask,(x0,x1))

        plt.imshow(s_mask,interpolation=None,cmap='jet')
        plt.show()
        print s_mask.shape
        exit()

        #Old method, sky_compress
        #if no mask
        #sky = np.median(image.data.reshape(ny,blocksize,nx,blocksize)
        #                .swapaxes(1,2).reshape(ny,nx,blocksize*blocksize),
        #                axis=2)
        #if mask
        #    data = np.ma.array(image.data, mask= (image.mask & bitmask),
        #                        fill_value=-1.)
        #    sky = np.ma.median(data.reshape(ny,blocksize,nx,blocksize)\
        #                       .swapaxes(1,2)\
        #                       .reshape(ny,nx,blocksize*blocksize), axis=2)
        #    sky = np.ma.getdata(sky)

    def ccd_mask(self,header_again=False,edge=15):
        '''Method to construct a mask, using the dimensions and positions
        of the CCDs
        Inputs:
        - again: whether to make do again the FITS header compilation of
        CCD dimensions/borders
        - edge: number of pixels to be removed from the edge of each CCD
        '''
        full_dim = CCD_position(rerun=header_again).tot_dim
        ccd_dim = CCD_position(rerun=header_again).ccd_dim
        print ccd_dim
        #check for higher/lower coordinates without data outside them
        for x in ccd_dim:
            if (x[1] > x[2]) or (x[3] > x[4]):
                logging.warning("Indices: initial > final")
                exit(1)
        n1 = (np.min(ccd_dim[:,1]),np.max(ccd_dim[:,3]))
        n2 = (np.min(ccd_dim[:,2]),np.max(ccd_dim[:,4]))
        #create array for FP without unused borders
        d1,d2 = np.ptp(n1),np.ptp(n2)
        arr = np.zeros((d1,d2))
        for ccd in ccd_dim:
            arr[ccd[1]-min(n1)-1+edge:ccd[2]-min(n1)-edge,
                ccd[3]-min(n2)-1+edge:ccd[4]-min(n2)-edge] = 1
        return arr

    def scale_mask(self,base,target_dim):
        '''Method adapted from: scalMask_dflat.py
        Receives a base mask and emulates its shape on the
        target array, using geometric scaling and nearest neighbor
        interpolation
        Inputs:
        - base: array to be redimensioned
        - target_dim: tuple with the dimensions of the target array
        '''
        #coordinates of the base mask
        x,y = np.arange(0,base.shape[1]),np.arange(0,base.shape[0])
        #defines an interpolator object
        intObj = scipy.interpolate.RegularGridInterpolator((y,x),
                                                        base.astype(int))
        #principle: identity matrix scaling with different scale for x and y
        #It's important to use the value of shape-1, otherwise the maximum
        #values are outside the base/target scaling dimensions
        ratio0 = (base.shape[0]-1.)/np.float(target_dim[0]-1.)
        ratio1 = (base.shape[1]-1.)/np.float(target_dim[1]-1.)
        #by multiplying target arrays indices by this ratios
        #will get the equivalent in the base mask dimension
        xaux,yaux = np.arange(0,target_dim[1]),np.arange(0,target_dim[0])
        xaux,yaux = ratio1*xaux,ratio0*yaux
        #at this point, just as double check, see if the maximum value is
        #out of bounds
        mssg = "Scaling mask: dimensions don't match"
        if (max(xaux) > base.shape[1]-1):
            logging.warning(mssg+" (dim 1:{0}, {1})".format(np.max(xaux),
                        base.shape[1]-1))
            xaux[-1] = (base.shape[1]-1.)
        if (max(yaux) > base.shape[0]-1):
            logging.warning(mssg+" (dim 0:{0}, {1})".format(np.max(yaux),
                        base.shape[0]-1))
            yaux[-1] = (base.shape[0]-1.)
        #construct the grid for interpolation
        xaux,yaux = np.meshgrid(xaux,yaux)
        cooaux = np.array(zip(yaux.ravel(),xaux.ravel()))
        #use nearest neighbor to find values
        fxy = intObj(cooaux,method='nearest').astype(bool)
        Maux = fxy.reshape(target_dim)
        #clean the interpolator
        intObj = None
        return Maux



if __name__=='__main__':
    #pos = CCD_position()
    Mimic().join_binned()
