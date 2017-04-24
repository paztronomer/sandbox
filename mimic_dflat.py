'''
Script to create mimic dome flats, by seeting up a focal plane array with
some feature and then applying the mask of the CCD sections, which effectively
catches the light
Different to my usual codes, I'll use self instead of cls 
'''
import os
import sys
import time
import numpy as np
#from scipy import ndimage
import tables
import fitsio
import pickle
#setup for display visualization
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

'''
STEPS
1) read and pickle header from an exposure, to get the CCDs positions --DONE
2) mask the borders of each CCD
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

class Mimic(CCD_position):
    def feature(self):
        back = None
        return back

    def mask(self):
        '''Method to construct a mask, using the dimensions and positions 
        of the CCDs 
        '''
        again = False
        edge = 15
        d1,d2 = CCD_position(rerun=again).tot_dim
        arr = np.zeros((d1,d2))
        for ccd in CCD_position(rerun=again).ccd_dim:
            arr[ccd[1]-1:ccd[2]+edge,ccd[3]-1:ccd[4]+edge] = 1

if __name__=='__main__':
    #pos = CCD_position()

    Mimic().mask()

