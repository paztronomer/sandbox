'''
Script to create mimic dome flats, by seeting up a focal plane array with
some feature and then applying the mask of the CCD sections, which effectively
catches the light
Different to my usual codes, I'll use method bounded to instances (self),
instead of methods bounded to classes (cls)
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
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

'''
STEPS
1) read and pickle header from an exposure, to get the CCDs positions --DONE
2) mask the borders of each CCD --DONE
3) create features --DONE
4) resize Fp and mask to different scales --DONE
Pending:
1) crop the borders, np.all(a == a[0,:], axis = 0) --DONE
'''

class Toolbox():
    def gauss2d(self,x1,x2,mu_x1,mu_x2,sigma_x1,sigma_x2):
        '''Method returns a 2D Gaussian, based in given vectors, means
        and standard deviations.
        Inputs:
        - x1,x2: 1D arrays to be used to define the meshgrid for which the
        Gaussian will be calculated.
        - mu_x1, mu_x2: float, means of both axis for which the Gaussian is
        defined
        - sigma_x1,sigma_x2: standard deviations for both the axis
        Outputs:
        - 2D array with the Gaussian profile
        '''
        x1,x2 = np.meshgrid(x1,x2)
        f1 = np.exp(-0.5*np.power((x1-mu_x1)/sigma_x1,2))
        f2 = np.exp(-0.5*np.power((x2-mu_x2)/sigma_x2,2))
        f = f1*f2/(2.*np.pi*sigma_x1*sigma_x2)
        return f

    def range_mean(self,arr,percent_low,percent_up):
        '''Method to calculate the mean of a range of values, given a
        distribution.
        Inputs:
        - arr: ndarray for which the mean will be calculated
        - percent_low,percent_up: limits of the range for which the mean will
        be calculated. Values in percent, from 0 t0 100
        Outputs:
        - value
        '''
        arr = arr.ravel()
        v1,v2 = np.percentile(arr,percent_low),np.percentile(arr,percent_up)
        tmp = arr[np.where(np.logical_and(arr>=v1,arr<=v2))]
        return np.mean(tmp)


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
        Inputs:
        - path: directory of the raw exposure (can be anyone)
        - fname: filename raw exposure (can be compressed or uncompressed)
        - rerun: wether to create again the array with the ccd limits and
        focal plane dimension
        Outputs:
        - assigns ccd limits to self.ccd_dim, and total focal plane dimensions
        to self.tot_dim. To do that it saves or loads pickle objects
        (ccd_dim.pickle,fp_dim.pickle)
        '''
        o1,o2 = 'ccd_dim.pickle','fp_dim.pickle'
        if rerun:
            d = os.path.join(path,fname)
            M_hdu = fitsio.FITS(d)
            dim = np.empty((0,5))
            for M in M_hdu[1:]:
                x = M.read_header()
                if ((x['CCDNUM'] <= 62) and (x['CCDNUM'] > 0)):
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
    def ccd_mask(self,header_again=False,edge=15,mask_again=False,
                exclude=[31]):
        '''Method to construct a mask, using the dimensions and positions
        of the CCDs. Converts it to sparse matrix, to save to disk: this
        is about 100 times lighter than normal numpy arrays (215 Mb vs 20 Gb)
        Inputs:
        - header_again: whether to make do again the FITS header compilation of
        CCD dimensions/borders
        - edge: number of pixels to be removed from the edge of each CCD
        - mask_again: wether to use the header info stored in the NPZ sparse
        matrix file, to re-construct the FP mask
        - exclude: list containig science ccdnum to be excluded from the mask
        Outputs:
        - sparse BSR matrix converted to dense array
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
                if ccd[0] not in exclude:
                    arr[ccd[1]-min(n1)-1+edge:ccd[2]-min(n1)-edge,
                        ccd[3]-min(n2)-1+edge:ccd[4]-min(n2)-edge] = 0
            #convert to sparse matrix
            bsr_arr = sparse.bsr_matrix(arr)
            sparse.save_npz(o1,bsr_arr)
        else:
            #load bsr matrix
            bsr_arr = sparse.load_npz(o1)
        #transform to dense array
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
        Outputs:
        - redimensioned array
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

    def get_mask(self,binsize,mask_again=False,header_again=False,
                exclude=[31]):
        '''Method to return the mask for the selected binning. A copy of
        the mask at resolution 16x16 is stored to speed up the loading.
        If other binsize is given, the mask will be recalculated.
        If the parameter mask_again is True, the mask will be computed from
        the header information
        The mask has 1 for outside region and 0 for inside CCDs
        Inputs:
        - binsize: tuple of 2 elements, giving the shape of the used bin
        - mask_again: wether to recalculate the mask from the header
        information, stored in ccd_dim.pickle, regardless the same binning
        may be used. If False, a stored mask of binning 16x16 will be loaded.
        - header_again: wether to recalculate the mask, reading again the
        header information and re-writting ccd_dim.pickle.
        - exclude: list of ccdnums to be excluded in the call of ccd_mask()
        method
        Outputs:
        - 2D array of 1/0 being the mask, which size defined by the binning
        '''
        if (binsize == (16,16) and not mask_again):
            bin_mask = np.load('s_mask_b{0}{1}.npy'.format(*binsize))
        elif (binsize != (16,16) and not mask_again):
            mask = Mimic().ccd_mask(mask_again=mask_again,
                                header_again=header_again)
            x0 = np.round(mask.shape[0]/np.float(binsize[0])).astype(int)
            x1 = np.round(mask.shape[1]/np.float(binsize[1])).astype(int)
            bin_mask = Mimic().scale_mask(mask,(x0,x1))
            np.save('s_mask_b{0}{1}.npy'.format(*binsize),bin_mask)
        elif mask_again:
            mask = Mimic().ccd_mask(mask_again=mask_again,
                                header_again=header_again)
            x0 = np.round(mask.shape[0]/np.float(binsize[0])).astype(int)
            x1 = np.round(mask.shape[1]/np.float(binsize[1])).astype(int)
            bin_mask = Mimic().scale_mask(mask,(x0,x1))
            np.save('s_mask_b{0}{1}.npy'.format(*binsize),bin_mask)
        else:
            logging.warning('Error loading the mask')
            return None
        return bin_mask

    def feature_multivariate(self,dim):
        '''Method to construct a test multivariate random feature. The seed is
        given so the result will be consistent between runs.
        I tuned the parameters to give an illumination pattern with a peak in a
        corner and two lower peaks along an edge
        After the multivariate sampling, a Gaussina filter is applied. Is the
        combination of this two methods who give the result.
        Inputs:
        - dim: fimensions of the binned focal plane
        Outputs:
        - a 2D array of the feature, with lower value being zero
        '''
        data = np.ones(dim)
        x0,x1 = dim
        mean = [x0/4,x1/3]
        cov = [[x0/16,0],[0,x1/2]]
        np.random.seed(20170425)
        #check the outputs from multivariate_normal, not sure about the meaning
        mvar = np.random.multivariate_normal(mean,cov,data.shape)
        data += mvar[:,:,0]
        data += mvar[:,:,1]
        aux = scipy.ndimage.filters.gaussian_filter(
            data,sigma=(x0/10,x1/3),order=(0,3),mode='reflect')
        aux += np.abs(np.min(aux))
        avg = Toolbox().range_mean(aux,90,100)
        aux = aux*(0.001/avg)
        return aux

    def feature1(self,dim):
        '''Feature-1, 2D gaussian, large pattern, vertical on the edge.
        Parameters were tunned by hand.
        Remember FWHM ~ 2.3548 sigma
        Inputs:
        - dim: dimensions of the binned focal plane
        Outputs:
        - 2D array with the feature
        '''
        y1,y2 = np.arange(0,dim[1]),np.arange(0,dim[0])
        #centers for the Gaussians
        mean_y1,mean_y2 = 0.5*dim[1],1.05*dim[0]
        #standard deviations for the distribution
        sigma_y1,sigma_y2 = dim[1]/2.35,0.8*dim[0]/2.35
        data = Toolbox().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        avg = Toolbox().range_mean(data,90,100)
        data = data*(0.001/avg)
        return data

    def feature2(self,dim):
        '''Feature-1, 2D gaussian, large pattern, rotated in 45 deg, maximum
        at the corner. To create this feature, a bigger array was defined, a
        Gaussian was included and then rotated. The returned array is a
        subsample of the original rotated, picked from the lower center.
        Parameters were tunned by hand.
        Remember FWHM ~ 2.3548 sigma
        Inputs:
        - dim: dimensions of the binned FP
        Outputs:
        - 2D array with the feature
        '''
        d0,d1 = 3*dim[0],3*dim[1]
        y1,y2 = np.arange(0,d1),np.arange(0,d0)
        #centers for the Gaussians
        mean_y1,mean_y2 = d1/5.,d0/3.
        #standard deviations for the distribution
        sigma_y1,sigma_y2 = .4*d1/2.35,2.*d0/2.35
        data = Toolbox().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        #rotate
        r_data = scipy.ndimage.interpolation.rotate(data,angle=-45.,axes=(1,0),
                                                    reshape=False,order=3,
                                                    mode='constant')
        avg = Toolbox().range_mean(r_data,90,100)
        r_data = r_data*(0.001/avg)
        return r_data[d0/3:d0/3+dim[0],d1/3:d1/3+dim[1]]

    def feature3(self,dim):
        '''Feature-1, 2D gaussian, 4 spots at the corners. The feature was
        constructed by adding Gaussians at the corners, one by one. Then, the
        last Gaussian probably has higher values, given the residuals of the
        previous ones. Parameters were tunned by hand.
        Remember FWHM ~ 2.3548 sigma
        Inputs:
        - dim: dimensions of the binned focal plane
        Output:
        - 2D array with the multiple features
        '''
        y1,y2 = np.arange(0,dim[1]),np.arange(0,dim[0])
        sigma_y1,sigma_y2 = 0.4*dim[1]/2.35,0.4*dim[0]/2.35
        #corner 1
        mean_y1,mean_y2 = 0*dim[1],0*dim[0]
        data = Toolbox().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        #corner 2
        mean_y1,mean_y2 = 0*dim[1],dim[0]
        data += Toolbox().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        #corner 3
        mean_y1,mean_y2 = dim[1],0*dim[0]
        data += Toolbox().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        #corner 4
        mean_y1,mean_y2 = dim[1],dim[0]
        data += Toolbox().gauss2d(y1,y2,mean_y1,mean_y2,sigma_y1,sigma_y2)
        avg = Toolbox().range_mean(data,90,100)
        data = data*(0.001/avg)
        return data

    def join_binned(self,binsize=16,header_again=False,mask_again=False,
                exclude=[31]):
        """Method to put together mask and data. To not mess with combining
        data from mask and image (averaging borders with data), both of them
        will be processed separately and then joined. Each feature will be
        saved in a separate numpy file
        Inputs:
        - binsize: if is an integer, use square bin, if a tuple, use its
        dimension
        - mask_again: wether to recalculate the mask from the header
        information, stored in ccd_dim.pickle, regardless the same binning
        may be used. If False, a stored mask of binning 16x16 will be loaded.
        - header_again: wether to recalculate the mask, reading again the
        header information and re-writting ccd_dim.pickle.
        - exclude: list of ccdnums to be excluded in the call of ccd_mask()
        Outputs:
        - none, it saves on files the features, with the applied mask
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
        s_mask = Mimic().get_mask(binsize,header_again=header_again,
                                mask_again=mask_again)
        #the feature data
        #the input focal plane image will be already in the new dimensions,
        #to save computational time
        mvar = Mimic().feature_multivariate(s_mask.shape)
        aux_minvar = np.min(mvar)
        mvar[s_mask.astype(bool)] = -1
        #
        feat1 = Mimic().feature1(s_mask.shape)
        aux_min1 = np.min(feat1)
        feat1[s_mask.astype(bool)] = -1
        #
        feat2 = Mimic().feature2(s_mask.shape)
        aux_min2 = np.min(feat2)
        feat2[s_mask.astype(bool)] = -1
        #
        feat3 = Mimic().feature3(s_mask.shape)
        aux_min3 = np.min(feat3)
        feat3[s_mask.astype(bool)] = -1
        if True:
            np.save('featmvar_b{0}{1}.npy'.format(*binsize),mvar)
            np.save('feat1_b{0}{1}.npy'.format(*binsize),mvar)
            np.save('feat2_b{0}{1}.npy'.format(*binsize),mvar)
            np.save('feat3_b{0}{1}.npy'.format(*binsize),mvar)
        if True:
            #checking by plot
            im = plt.imshow(feat3,cmap='viridis',interpolation='none',
                        origin='upper',vmin=aux_min3)
            plt.colorbar(im)
            plt.show()
        print '\nFeatures for binning {0}x{1} were saved'.format(*binsize)
        return True


if __name__=='__main__':
    #pos = CCD_position()
    Mimic().join_binned()
