"""Simple script to count exposures for bias and flats, for a specific nite
"""

import os
import time
import subprocess
import pwd
import numpy as np
import despydb.desdbi as desdbi
import fitsio 
import pickle

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
    def calib_nite(cls,nite):
        q = "select expnum,band,obstype,nite,time_obs from exposure"
        q += " where nite={0}".format(nite)
        q += " and obstype in ('dome flat','zero')"
        q += " order by expnum"
        datatype = ['i4','a10','a50','i4','a50']
        return Toolbox.dbquery(q,datatype)

class Count():
    @classmethod
    def count(cls,nite,exclude=None):
        tab = Toolbox.calib_nite(nite)
        #here there is no need for remove expnum duplicates
        s1 = []
        for lab in np.unique(tab['obstype']):
            for b in np.unique(tab['band']):
                a1 = tab[(tab['obstype']==lab) & (tab['band']==b)].shape[0]
                if exclude:
                    tab_aux = tab[~np.in1d(tab['expnum'],exclude)]
                    a2 = tab_aux[(tab_aux['obstype']==lab) & 
                                (tab_aux['band']==b)].shape[0]
                else:
                    a2 = np.nan
                s1.append([lab,b,a1,a2])
        print '================================================='
        print '{0:<12}{1:<12}{2:<8}'.format(*['OBSTYPE','BAND','N_EXP'])
        for i in xrange(len(s1)):
            print '{0:<12}{1:<12}{2:<8}{3:<8}'.format(*s1[i])
        print '================================================='
        return s1

if __name__=='__main__':
    print time.ctime()
    Count.count(20170206,exclude='y4_e2.BAD')
