"""Simple script with a couple of tools for structured arrays
"""
import os
import sys
import time 
import numpy as np
import pandas as pd
import despydb.desdbi as desdbi

#########################DES##############################
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
    def aux(cls):
        q = "select nite,band,expnum"
        q += " from exposure"
        q += " where nite between 20170201 and 20170223"
        q += " and obstype='dome flat'"
        q += " and band='r'"
        datatype = ['i4','a10','i4']
        tab = Toolbox.dbquery(q,datatype)
        return tab
##########################################################

class Use():
    @classmethod
    def sarr_sort(cls,arr,column):
        '''Inputs
        - arr: structured array to be sorted
        - column: column to sort
        Returns
        - sorted structured array
        '''
        return np.sort(arr,order=column)

    @classmethod
    def sarr_drop(cls,arr,column,posit='first'):
        """Method (based on np.unique) to remove duplicates, keeping first
        or last row. Recommend to sort the structured array before the
        drop.
        Inputs:
        - arr: structured array
        - column: column to be used as drop criteria
        - posit: accepts either first or last, to keep the 1st/last occurence
        Outputs:
        - structured array 
        """
        #option1
        #np.unique returns: unique array,indices of the original array to be 
        #used to final array,indices of the unique array that can be used to 
        #reconstruct the original,number of times each of the unique array 
        #elements appear in the original
        posit = posit.lower()
        uarr,uidx,uinverse,ucounts= np.unique(
            arr[column],return_index=True,
            return_inverse=True,
            return_counts=True)
        if posit == 'first':
            #by default np.unique picks the first
            out = arr[uidx]
        elif posit == 'last':
            lastidx = []
            #using the amount of points of the unique array, pick the last 
            #occurrence in the array of indices
            for idx in xrange(uarr.shape[0]):
                lastidx.append(np.where(uinverse==idx)[0][-1])
            out = arr[np.array(lastidx)]
        else:
            out = False
            m = 'Please give a valid position to be selected.'
            m += 'Valid values are: first, last'
            raise ValueError(m)
        returns out

if __name__ == '__main__':
    
    #tab = Toolbox.aux()
    #np.save('tab.npy',tab)
    tab = np.load('tab.npy')
    Use.sarr_drop(tab,'nite')
