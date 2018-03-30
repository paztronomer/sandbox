'''Script to create a average-like image from a set, using the median
'''

import os
import time
import numpy as np
import fitsio
import pickle
import logging
import argparse
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
#
#

def db_red_pixcor(arr_expnum, reqnum, band,
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

def main_aux(explist=None, pathlist=None, reqnum=None, band=None):
    if ((explist is not None) and (pathlist is None)):
        # Load the table
        exp = np.genfromtxt(
            explist,
            dtype={'names' : ['expnum'], 'formats' : ['i8']},
            comments='#',
            missing_values=np.nan,
            usecols=0
        )
        # Get the full path from the DB, get in chunks of 1000 (max before
        # use temp tables)
        if (exp.size > 999):
            accum = []
            logging.warning('Still not implemented!')
            exit(1)
        fpath = db_red_pixcor(exp['expnum'], reqnum, band)
    elif ((pathlist is not None) and (explist is None)):
        # Load the table, with no constraint on filetype, in case path
        # is extremely large
        exp = np.genfromtxt(
            explist,
            names = ['path'],
            comments='#',
            missing_values=np.nan,
            usecols=0
        )
        fpath = exp['path']
    # Remove duplicates
    uarr, uidx, uinverse, ucounts= np.unique(
        fpath,
        return_index=True,
        return_inverse=True,
        return_counts=True,
    )
    # With the full paths open the files in parallel, read in boxes of
    # 256 sq pix, stacking them
    h = fpath[0]
    x0, y0, x1, y1 = 0, 0, 256, 256
    s = fits_section(h, x0, y0, x1, y1)
    print(s)
    # Go pixel by pixel doing the median

def fits_section(fname, x0, y0, x1, y1, ext=0, bin=256):
    ''' Load fits file section bin size of sq pixels
    '''
    if (int(abs(x1 - x0)) == bin) and (int(abs(y1 - y0)) == bin):
        pass
    else:
        logging.warning('Coodinates don\'t have the size of the input bin')
    if os.path.exists(fname):
        fits = fitsio.FITS(fname)
        sq = np.copy(fits[ext][y0:y1 , x0:x1])
        fits.close()
        return sq
    else:
        logging.error('File {0} does not exists'.format(fname))
        return None

if __name__ == '__main__':
    txt_gral = 'Code to create a median image from an input list of expnum'
    txt_gral += ' or full paths. No normalization is applied so make sure'
    txt_gral += ' all images are on the same ground level'
    par = argparse.ArgumentParser(description=txt_gral)
    src = par.add_mutually_exclusive_group()
    #
    exp_tmp = '/Users/fco/Code/des_calibrations/skytemplates_build'
    exp_tmp += '/Y5_skytemplates/expnum_g_not20171109t1125.csv'
    h0 = ''
    src.add_argument('--explist', help=h0, default=exp_tmp)
    h1 = ''
    src.add_argument('--pathlist', help=h1)
    h2 = 'Reqnum in case --explist was feeded'
    par.add_argument('--req', help=h2, type=int)
    h3 = 'Band in case --explist was feeded'
    par.add_argument('--band', help=h3, type=str)
    # Parse args
    par = par.parse_args()
    in_kw = {
        'explist' : par.explist,
        'pathlist' : par.pathlist,
        'reqnum' : par.req,
        'band' : par.band,
    }
    main_aux(**in_kw)
