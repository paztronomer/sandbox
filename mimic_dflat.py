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
import scipy.ndimage
import scipy.sparse as sparse
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
    def gauss2d(self,x1,x2,mu_x1,mu_x2,sigma_x1,sigma_x2):
        '''Method returns a 2D Gaussian, based in different means
        and standard deviations
        '''
        x1,x2 = np.meshgrid(x1,x2)
        f1 = np.exp(-0.5*np.power((x1-mu_x1)/sigma_x1,2))
        f2 = np.exp(-0.5*np.power((x2-mu_x2)/sigma_x2,2))
        f = f1*f2/(2.*np.pi*sigma_x1*sigma_x2)
        f = f.T
        return f

    def feature_multivariate(self,dim):
        '''Method to construct a test multivariate feature. The seed is given
        so the result will be consistent between runs.
        I tuned the parameters to give an illumination pattern with a peak in a
        corner and two lower peaks along an edge
        After the multivariate sampling, a Gaussina filter is applied. Is the
        combination of this two methods who give the results and not alone.
        Inputs:
        - dim: fimensions of the binned FP
        Returns:
        - a 2D array with lower value being zero
        '''
        data = np.ones(dim)
        x0,x1 = dim
        mean = [x0/4,x1/3]
        cov = [[x0/16,0],[0,x1/2]]
        np.random.seed(20170425)
        mvar = np.random.multivariate_normal(mean,cov,data.shape)
        data += mvar[:,:,0]
        data += mvar[:,:,1]
        aux = scipy.ndimage.filters.gaussian_filter(
            data,sigma=(x0/10,x1/3),order=(0,3),mode='reflect')
        aux += np.abs(np.min(aux))
        return aux

    def feature1(self,dim):
        '''Feature-1, 2D gaussian, large pattern, vertical on the edge.
        Remember FWHM ~ 2.3548 sigma
        Inputs:
        - dim: dimensions of the binned FP
        Output:
        - 2D array with the feature
        '''
        y1,y2 = np.arange(0,dim[0]),np.arange(0,dim[1])
        #centers for the Gaussians
        mean_y1,mean_y2 = 0.5*dim[0],1.05*dim[1]
        #standard deviations for the distribution
        sigma_y1,sigma_y2 = dim[0]/2.35,0.8*dim[1]/2.35
        data = Mimic().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        return data

    def feature2(self,dim):
        '''Feature-1, 2D gaussian, large pattern, rotated in 45 deg, maximum
        at the corner.
        Remember FWHM ~ 2.3548 sigma
        Inputs:
        - dim: dimensions of the binned FP
        Output:
        - 2D array with the feature
        '''
        d0,d1 = 3*dim[0],3*dim[1]
        y1,y2 = np.arange(0,d0),np.arange(0,d1)
        #centers for the Gaussians
        mean_y1,mean_y2 = d0/5.,d1/3.
        #standard deviations for the distribution
        sigma_y1,sigma_y2 = .4*d0/2.35,2.*d1/2.35
        data = Mimic().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        #rotate
        r_data = scipy.ndimage.interpolation.rotate(data,angle=-45.,axes=(1,0),
                                                    reshape=False,order=3,
                                                    mode='constant')
        return r_data[d0/3:d0/3+dim[0],d1/3:d1/3+dim[1]]

    def feature3(self,dim):
        '''Feature-1, 2D gaussian, 4 spots at the corners.
        Remember FWHM ~ 2.3548 sigma
        Inputs:
        - dim: dimensions of the binned FP
        Output:
        - 2D array with the multiple features
        '''
        y1,y2 = np.arange(0,dim[0]),np.arange(0,dim[1])
        #centers for the Gaussians
        mean_y1,mean_y2 = 0*dim[0],0*dim[1]
        #standard deviations for the distribution
        sigma_y1,sigma_y2 = 0.4*dim[0]/2.35,0.4*dim[1]/2.35
        data = Mimic().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        #
        mean_y1,mean_y2 = 0*dim[0],dim[1]
        data += Mimic().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        #
        mean_y1,mean_y2 = dim[0],0*dim[1]
        data += Mimic().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        #
        mean_y1,mean_y2 = dim[0],dim[1]
        data += Mimic().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        return data

    def get_mask(self,binsize,mask_again=False):
        '''Method to return the mask for the selected binning. A copy of
        the mask at resolution 16x16 is stored to speed up the loading.
        If other binsize is given, the mask will be recalculated.
        If the parameter mask_again is True, the mask will be computed from
        the header information
        The mask has 1 for outside region and 0 for inside CCDs
        Inputs:
        - binsize: tuple of 2 elements, giving the shape of the used bin
        - mask_again: wether to recalculate the mask from the header
        information, stored in ccd_dim.pickle. This will call ccd_mask()
        method
        Returns:
        - mask in the size defined by the binning
        '''
        if (binsize == (16,16) and not mask_again):
            bin_mask = np.load('s_mask_b16.npy')
        elif (binsize != (16,16) and not mask_again):
            mask = Mimic().ccd_mask(mask_again=mask_again)
            x0 = np.round(mask.shape[0]/np.float(binsize[0])).astype(int)
            x1 = np.round(mask.shape[1]/np.float(binsize[1])).astype(int)
            bin_mask = Mimic().scale_mask(mask,(x0,x1))
        elif mask_again:
            mask = Mimic().ccd_mask(mask_again=mask_again)
            x0 = np.round(mask.shape[0]/np.float(binsize[0])).astype(int)
            x1 = np.round(mask.shape[1]/np.float(binsize[1])).astype(int)
            bin_mask = Mimic().scale_mask(mask,(x0,x1))
        else:
            logging.warning('Error loading the mask')
            return None
        return bin_mask

    def join_binned(self,binsize=16,header_again=False):
        """Method to put together mask and data. To not mess with combining
        data from mask and image (averaging borders with data), both of them
        will be processed separately and then joined
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

        #the mask
        #as mask is a binary-valued object, nearest neighbor interpolation is
        #the selected algorithm
        s_mask = Mimic().get_mask(binsize)

        #the feature data
        #the input focal plane image will be already in the new dimensions,
        #to save computational time
        if False:
            mvar = Mimic().feature_multivariate(s_mask.shape)
            aux_min = np.min(mvar)
            mvar[s_mask.astype(bool)] = -1
        if False:
            feat1 = Mimic().feature1(s_mask.shape)
            aux_min = np.min(feat1)
            feat1[s_mask.astype(bool)] = -1
        if False:
            feat2 = Mimic().feature2(s_mask.shape)
            aux_min = np.min(feat2)
            feat2[s_mask.astype(bool)] = -1
        if False:
            feat3 = Mimic().feature3(s_mask.shape)
            aux_min = np.min(feat3)
            feat3[s_mask.astype(bool)] = -1
            #checking by plot
            im = plt.imshow(feat3,cmap='viridis',interpolation=None,
                        origin='lower',vmin=aux_min)
            plt.colorbar(im)
            plt.show()
        return True

    def ccd_mask(self,header_again=False,edge=15,mask_again=True):
        '''Method to construct a mask, using the dimensions and positions
        of the CCDs. Converts it to sparse matrix, to save to disk: this
        is about 100 times lighter than normal numpy arrays (215 Mb vs 20 Gb)
        Inputs:
        - again: whether to make do again the FITS header compilation of
        CCD dimensions/borders
        - edge: number of pixels to be removed from the edge of each CCD
        - mask_again: wether to use the header info again to construct the
        FP mask again
        '''
        o1 = "ccd_mask_bsr.npz"
        if mask_again:
            full_dim = CCD_position(rerun=header_again).tot_dim
            ccd_dim = CCD_position(rerun=header_again).ccd_dim
            #check for higher/lower bsrrdinates without data outside them
            for x in ccd_dim:
                if (x[1] > x[2]) or (x[3] > x[4]):
                    logging.warning("Indices: initial > final")
                    exit(1)
            n1 = (np.min(ccd_dim[:,1]),np.max(ccd_dim[:,2]))
            n2 = (np.min(ccd_dim[:,3]),np.max(ccd_dim[:,4]))
            #create array for FP without unused borders
            d1,d2 = np.ptp(n1),np.ptp(n2)
            arr = np.ones((d1,d2))
            for ccd in ccd_dim:
                arr[ccd[1]-min(n1)-1+edge:ccd[2]-min(n1)-edge,
                    ccd[3]-min(n2)-1+edge:ccd[4]-min(n2)-edge] = 0
            #convert to sparse matrix
            bsr_arr = sparse.bsr_matrix(arr)
            sparse.save_npz(o1,bsr_arr)
        else:
            #load bsr matrix
            bsr_arr = sparse.load_npz(o1)
            if not isinstance(bsr_arr,np.ndarray):
                bsr_arr = bsr_arr.toarray()
        return bsr_arr

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
