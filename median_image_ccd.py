'''Script to create a average-like image from a set, using the median. No 
class-structure for an easy call in parallel
'''

import os
import time
import uuid
import numpy as np
import pandas as pd
import fitsio
import pickle
import logging
import argparse
try:
    import partial
except:
    logging.warning('No available module: partial')
import multiprocessing as mp
import easyaccess as ea

# Setup logging
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO,
)



class Toolbox():
    @classmethod
    def dbquery(cls,toquery,outdtype,dbsection='db-desoper',help_txt=False):
        '''the personal setup file .desservices.ini must be pointed by desfile
        DB section by default will be desoper
        '''
        desfile = os.path.join(os.getenv('HOME'),'.desservices.ini')
        section = dbsection
        dbi = desdbi.DesDbi(desfile,section)
        if help_txt: help(dbi)
        cursor = dbi.cursor()
        cursor.execute(toquery)
        cols = [line[0].lower() for line in cursor.description]
        rows = cursor.fetchall()
        outtab = np.rec.array(rows,dtype=zip(cols,outdtype))
        return outtab

    @classmethod
    def niterange_supercal(cls,niterange,reqnum):
        N1,N2 = niterange
        query = 'select i.path, n.filename, n.band, n.nite, n.expnum,'
        query += ' att.reqnum,att.attnum'
        query += ' from flat_qa q, miscfile m, miscfile n,file_archive_info i,'
        query += ' pfw_attempt att'
        query += ' where m.nite between {0} and {1}'.format(N1,N2)
        query += ' and att.reqnum={0}'.format(reqnum)
        query += ' and m.pfw_attempt_id=n.pfw_attempt_id'
        query += ' and att.id=n.pfw_attempt_id'
        query += ' and m.expnum=n.expnum'
        query += ' and q.filename=m.filename'
        query += ' and n.filename=i.filename'
        datatype = ['a100','a100','a5','i4','i4','i4','i4']
        tab = Toolbox.dbquery(query,datatype)
        return tab


class FPBinned():
    def __init__(self,folder,fits):
        '''Simple method to open focal plane binned images
        When a position not belongs to focal plane, the value is -1
        Before return it, add 1 to set outer region to zero value
        '''
        fname = os.path.join(folder,fits)
        M_header = fitsio.read_header(fname)
        M_hdu = fitsio.FITS(fname)[0]
        tmp = M_hdu.read()
        #tmp = Toolbox.detect_outlier(tmp)
        tmp += 1.
        self.fpBinned = tmp


class Specific():
    @classmethod
    def weight_image(cls,db_table):
        '''Create weight images by the median of the distribution, pixel
        by pixel, and save them
        '''
        print( 'Columns of the info input table (from DESDB):'.format(
            db_table.dtype.names))
        #receives a DB table for a niterange, with excluded list already
        #applied, and  calculate the median image
        band = np.unique(db_table['band'])
        nite = np.unique(db_table['nite'])
        rootp = '/archive_data/desarchive'
        #select per nite per band
        for b in band:
            for idx,n in enumerate(nite):
                auxsel = db_table[(db_table['band']==b) &
                                (db_table['nite']==n)]
                for idx2,exp in enumerate(np.unique(auxsel['expnum'])):
                    '''#as multiple reqnum could be associted with a same exposure,
                    #we will select the highest reqnum to be used for each one
                    #of the exposures
                    reqmax = np.max(auxsel[(auxsel['expnum']==exp)]['reqnum'])
                    #now create the median of the night, per band
                    pth = auxsel[(auxsel['expnum']==exp) &
                                (auxsel['reqnum']==reqmax)]['path'][0]
                    fnm = auxsel[(auxsel['expnum']==exp) &
                                (auxsel['reqnum']==reqmax)]['filename'][0]
                    '''
                    #now create the median of the night, per band
                    pth = auxsel[auxsel['expnum']==exp]['path'][0]
                    fnm = auxsel[auxsel['expnum']==exp]['filename'][0]
                    if idx2 == 0:
                        M = FPBinned(os.path.join(rootp,pth),fnm).fpBinned
                        aux_med = np.median(M)
                        #============================
                        #for compare_dflat_binned_fp
                        #must remove the median norm
                        #============================
                        M /=  aux_med
                    else:
                        N = FPBinned(os.path.join(rootp,pth),fnm).fpBinned
                        aux2_med = np.median(N)
                        #============================
                        #for compare_dflat_binned_fp
                        #must remove the median norm
                        #============================
                        N /= aux2_med
                        M = np.dstack((M,N))
            #call a method to get the median
            M_w = Specific.stat_cube(M,(lambda: np.median)())
            direc = 'weighted/'
            fnm = 'y4e2_{0}.fits'.format(b)
            fits = fitsio.FITS(os.path.join(direc,fnm),'rw')
            fits.write(M_w)
            fits[-1].write_checksum()
            fits.close()
        return True

    @classmethod
    def stat_cube(cls,arr3,f):
        '''Receives a data cube (3D array) and performs the given statistics
        over the third dimension, pixel by pixel
        Uses numpy iteration tools
        '''
        out = np.zeros_like(arr3[:,:,0])
        it = np.nditer(arr3[:,:,0], flags=['multi_index'])
        while not it.finished:
            i1,i2 = it.multi_index
            out[i1,i2] = f(arr3[i1,i2,:])
            it.iternext()
        return out

    @classmethod
    def norm_image(cls,db_table):
        '''Divides a set of images by the input median image.
        The weighted images must be already calculated. The single images to
        be normalized are assumed to be on disk, in the DES filesystem
        '''
        #receives a DB table for a niterange, calculate the median image
        band = np.unique(db_table['band'])
        nite = np.unique(db_table['nite'])
        rootp = '/archive_data/desarchive'
        name_list = []
        #select per nite per band per exposure
        for bb in band:
            for idx,nn in enumerate(nite):
                nlist = []
                auxsel = db_table[(db_table['band']==bb) &
                                (db_table['nite']==nn)]
                for idx2,exp in enumerate(np.unique(auxsel['expnum'])):
                    #with this, save filename/path to be used to feed a list
                    tmp = auxsel[auxsel['expnum']==exp][0]
                    nlist.append([tmp['path'],tmp['filename']])
                #the index of the exposure to be normalized and displayed
                #will increase for night to night and then stuck and
                #repeat
                if idx < len(nlist):
                    indsel = idx
                else:
                    indsel = len(nlist)-1
                #we have the index of the filename/path list, so now do the job
                #open the median-image
                norm = FPBinned('weighted','y4e2_{0}.fits'.format(bb)).fpBinned
                #open the selected pixcor image
                pixc = FPBinned(os.path.join(rootp,nlist[indsel][0]),
                            nlist[indsel][1]).fpBinned
                #============================
                #for compare_dflat_binned_fp
                #must remove the median norm
                #============================
                pixc /= np.median(pixc)
                #save the normalized exposure on ds9plot/
                res = pixc/norm
                direc = 'ds9plot_y4e2/'
                fnm = 'norm_{0}_{1}.fits'.format(bb,nn)
                fits2 = fitsio.FITS(os.path.join(direc,fnm),'rw')
                fits2.write(res)
                fits2[-1].write_checksum()
                fits2.close()
                #save the filename
                name_list.append(os.path.join(direc,fnm))
        return name_list


class Call:
    @classmethod
    def wrap1(cls,niterange,reqnum,exclusion=None):
        queryDB = False
        if queryDB:
            with open('niterange.pickle','wb') as aid:
                dbtab = Toolbox.niterange_supercal(niterange,reqnum)
                pickle.dump(dbtab,aid,protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open('niterange.pickle','rb') as aid:
                dbtab = pickle.load(aid)
        dbtab = dbtab[dbtab['nite']!=20170211]
        #remove the excluded list
        if exclusion:
            exclusion = np.loadtxt(exclusion)
            dbtab = dbtab[~np.in1d(dbtab['expnum'],exclusion)]
        #remove duplicates
        uarr,uidx,uinverse,ucounts= np.unique(
            dbtab['expnum'],return_index=True,
            return_inverse=True,
            return_counts=True)
        dbtab = dbtab[uidx]
        #create weighted image
        Specific.weight_image(dbtab)
        #normalize the images by its median
        Specific.norm_image(dbtab)
#
# ============================================================================
#
#

def db_red_pixcor(reqnum, band,
                  root_path='/archive_data/desarchive/'):
    ''' Simple function to get path information from DB, for a set under same
    reqnum and band
    '''
    connect = ea.connect('desoper')
    cursor = connect.cursor()
    q = 'select fai.path, fai.filename, fai.compression, e.band'
    q += ' from file_archive_info fai, desfile d, pfw_attempt att, exposure e'
    q += ' where fai.desfile_id=d.id'
    q += ' and d.filetype=\'red_pixcor\''
    q += ' and d.pfw_attempt_id=att.id'
    q += ' and att.reqnum={0}'.format(reqnum)
    q += ' and att.unitname=CONCAT(\'D00\', e.expnum)'
    q += ' and e.band=\'{0}\''.format(band)
    outtab = connect.query_to_pandas(q)
    # Transform column names to lower case
    outtab.columns = map(str.lower, outtab.columns)
    # Double check, remove duplicates
    outtab.drop_duplicates(['filename'], inplace=True)
    outtab.reset_index(drop=True, inplace=True)
    # Construct the full path
    def f1(a, b, c):
        if (c == None): return os.path.join(root_path, a, b)
        else: return os.path.join(root_path, a, b + c)
    # Construct the full path
    full_path = map(f1, outtab['path'], outtab['filename'],
                    outtab['compression'])
    full_path = list(full_path)
    full_path = np.array(full_path)
    return full_path

def fits_section(aux_list):
    ''' Load fits file section bin size of sq pixels
    '''
    fname, coo, ext, delta = aux_list
    x0, y0, x1, y1 = coo
    if (int(abs(x1 - x0)) == delta) and (int(abs(y1 - y0)) == delta):        
        pass
    else:
        logging.warning('Coordinates don\'t have the size of the input bin')
    if os.path.exists(fname):
        fits = fitsio.FITS(fname)
        sq = np.copy(fits[ext][y0:y1 , x0:x1])
        fits.close()
        return sq
    else:
        logging.error('File {0} does not exists'.format(fname))
        print('returning -1')
        return -1

def stat_cube(x3d, func):
    ''' Function to calculate the median image per pixel, using an input 
    3dimensional array containing the stamps on which to work.
    Uses numpy iteration tools
    '''
    out = np.zeros_like(x3d[:, :, 0])
    it = np.nditer(x3d[:, :, 0], flags=['multi_index'])
    while not it.finished:
        i1, i2 = it.multi_index
        out[i1, i2] = func(x3d[i1, i2, :])
        it.iternext()
    return out

def main_aux(pathlist=None, 
             reqnum=None, 
             band=None, 
             nproc=None, 
             chunk=None,
             px_side=None,
             label=None,):
    ''' Main auxiliaty function to call the code in parallel
    '''
    if (pathlist is None):
        # If no set of pull paths is provided, go to the DB and retrieve them,
        # based on reqnum and band
        fpath = db_red_pixcor(reqnum, band)
    elif (pathlist is not None):
        # Load the table
        exp = np.genfromtxt(
            pathlist,
            dtype=[('path','|S200'),],
            comments='#',
            missing_values=np.nan,
            usecols=0
        )
        fpath = exp['path']
    # Remove duplicates using np.unique capabilities
    uarr, uidx, uinverse, ucounts= np.unique(
        fpath,
        return_index=True,
        return_inverse=True,
        return_counts=True,
    )
    if False:
        # Save table for simplicity of testing
        tmp = pd.DataFrame({'path' : fpath})
        tmp.to_csv('path_r3437_g.txt', index=False, header=False)
    #
    # NOTE
    # The goal is to read the same section, from a bunch of images, in  
    # parallel. Then stack the set of section and calculate the median 
    # of each pixel.
    #
    # Get shape of the CCD
    try:
        aux_fits = fitsio.read(fpath[0])
        f_dim = aux_fits.shape
        logging.info('CCD shape: {0}'.format(f_dim))
    except:
        logging.error('File cannot be read. Using 2048x4096 as CCD dimensions')
        f_dim = (4096, 2048)
    # Define squared sections to be used for calculation
    if (px_side is None):
        px_side = 512
    if ((f_dim[0] % px_side != 0) or (f_dim[1] % px_side != 0)):
        logging.warning('{0}x{0} des not fit exactly on CCD'.format(px_side))
    # idx_d0 = np.arange(0, fits_dim[0] + 1, px_side)
    # idx_d1 = np.arange(0, fits_dim[1] + 1, px_side)
    yx = np.mgrid[0 : f_dim[0]+1 : px_side, 0 : f_dim[1]+1 : px_side]
    # coo_list contains [x0, y0, x1, y1]
    coo_list = []
    for iy in range(f_dim[0] / px_side):#yx.shape[1]):
        for ix in range(f_dim[1] / px_side):#(yx.shape[2]):
            y0, x0 = (yx[0, iy, ix], yx[1, iy, ix])
            y1, x1 = (yx[0, iy+1, ix], yx[1, iy, ix+1])
            coo_list.append((x0, y0, x1, y1))
    # Create an array to harbor the results
    tmp_res = np.zeros((f_dim[0], f_dim[1]))
    # Setup the parallel call
    # The value for chunk will help prevent memory errors. "Chops the iterable 
    # into a number of chunks which it submits to the process pool as separate 
    # tasks"
    if (nproc is None):
        nproc = mp.cpu_count()
    if (chunk is None):
        chunk = int( np.ceil(fpath.size / nproc) )
    logging.info('Launch {0} parallel processes'.format(nproc))
    P1 = mp.Pool(processes=nproc)
    # Here I can also parallelize in terms of coordinates, but need to be 
    # careful to not mix results from different sections of the CCD, and 
    # manage the Pool (2 in fact)
    for idx_c, coo in enumerate(coo_list):
        logging.info('Box {0} of {1}'.format(idx_c + 1, len(coo_list)))
        # coo: coordinates, ext: extension to read from the FITS file, 
        # delta: side size (pixel) of the square used to sample the FITS file 
        kw_section = {
            'coo' : coo,
            'ext' : 0,
            'delta' : px_side,
        }
        # Call using partial() or constructed list
        try:
            partial_aux = partial(
                fits_section, 
                [kw_section['coo'], kw_section['ext'], kw_section['delta']]
            )
            t0 = time.time()
            # map_async does not block the processes, and executes non-ordered
            boxi = P1.map_async(partial_aux, fpath, chunk)
            boxi.wait()  
            t1 = time.time()
        except:
            aux_list = [(fnm, 
                         kw_section['coo'], 
                         kw_section['ext'], 
                         kw_section['delta']) for fnm in fpath]
            t0 = time.time()
            boxi = P1.map_async(fits_section, aux_list, chunk)
            boxi.wait()
            t1 = time.time()
        boxi = boxi.get()
        # At this point the Pool has returned a list, containig the stamps
        # for all the CCDs
        # Call the median image generator, using np.nditer(), with a 3D array
        x3d = np.dstack(boxi)
        z_median = stat_cube(x3d, (lambda: np.median)())
        # Put into the auxiliary array 
        x0, y0, x1, y1 = coo
        tmp_res[y0:y1, x0:x1] = z_median
    # The results looks good!
    # Save in FITS format
    if (label is None):
        label = str(uuid.uuid4())
    ores = 'medImg_{0}.fits'.format(label)
    while os.path.exists(ores):
        logging.info('File {0} exists. Changing output name'.format(ores))
        ores = 'medImg_{0}.fits'.format(str(uuid.uuid4())) 
    fits = fitsio.FITS(ores,'rw')
    fits.write(tmp_res)
    fits[-1].write_checksum()
    fits.close()
    logging.info('Median image {0} saved'.format(ores))
    return True


if __name__ == '__main__':
    txt_gral = 'Code to create a median image from an input list of expnum'
    txt_gral += ' or full paths. No normalization is applied so make sure'
    txt_gral += ' all images are on the same ground level, example: red_pixcor'
    par = argparse.ArgumentParser(description=txt_gral)
    #
    t = 'path_r3437_g.txt'
    b = None #'g' 
    r = None #3437
    h1 = 'List of full path to the exposures. Default: {0}'.format(t)
    par.add_argument('--pathlist', help=h1, default=t)
    h2 = 'Reqnum in case no --pathlist was feeded. Default: {0}'.format(r)
    par.add_argument('--req', help=h2, type=int, default=r)
    h3 = 'Band in case no --pathlist was feeded. Default: {0}'.format(b)
    par.add_argument('--band', help=h3, type=str, default=b)
    h4 = 'Number or parallel processes to run. Default is number of CPUs'
    par.add_argument('-n', help=h4, metavar='', type=int)
    h5 = 'Chunks on which to divide the inputs in the processes. Default: 1'
    par.add_argument('-c', help=h5, metavar='', type=int)
    h6 = 'Label to use for output median image. Default is a random string'
    par.add_argument('--label', help=h6, metavar='')
    h7 = 'Side of the square (in pix) to be used as sections for dividing'
    h7 += ' the CCD. Need to fit exactly into the CCD array. Default: 512'
    par.add_argument('--square', help=h7, metavar='', type=int)
    # Parse args
    par = par.parse_args()
    in_kw = {
        'pathlist' : par.pathlist,
        'reqnum' : par.req,
        'band' : par.band,
        'nproc' : par.n,
        'chunk' : par.c,
        'label' : par.label,
        'px_side' : par.square,
    }
    main_aux(**in_kw)
