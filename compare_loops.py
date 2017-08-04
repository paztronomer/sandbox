"""Compare some loops
Francisco Paz-Chinchon
"""
import os
import sys
import time
import timeit
import gc
import numpy as np

class LC():
    @classmethod
    def open_lc(cls,fname):
        tmp = np.genfromtxt(fname,dtype=[np.float,np.float],comments="#",
                            filling_values=np.nan,usecols=(0,1),
                            names=["time","flux"])
        return tmp


class Play():
    """Class to experiment
    """
    def __init__(self,lc):
        self.lc = lc

    def rate1(self,xseed=20171224):
        data = self.lc
        np.random.seed(xseed)
        d0 = data.shape[0]
        pcent = 0.25
        rand_i = np.random.randint(1,high=d0-1,size=int(np.floor(d0*pcent)))
        sampl =[]
        for i in rand_i:
            sampl.append(data["time"][i]-data["time"][i-1])
            #compare with doing N += some, counter += 1
        sampl = np.mean(np.array(sampl))
        return sampl

    def rate2(self,xseed=20171224):
        data = self.lc
        np.random.seed(xseed)
        d0 = data.shape[0]
        pcent = 0.25
        rand_i = np.random.randint(1,high=d0-1,size=int(np.floor(d0*pcent)))
        S,c = 0,0
        for i in rand_i:
            S += data["time"][i]-data["time"][i-1]
            c += 1
        return S/c

    def rate3(self,xseed=20171224):
        """Optimizing loops
        https://www.python.org/doc/essays/list2str/ 
        """
        data = self.lc
        np.random.seed(xseed)
        d0 = data.shape[0]
        pcent = 0.25
        rand_i = np.random.randint(1,high=d0-1,size=int(np.floor(d0*pcent)))
        def D(j):
            return data["time"][j]-data["time"][j-1]
        aux_arr = map(D,rand_i)
        return np.sum(aux_arr)/len(aux_arr)

    def rate4(self,xseed=20171224):
        data = self.lc
        np.random.seed(xseed)
        d0 = data.shape[0]
        pcent = 0.25
        rand_i = np.random.randint(1,high=d0-1,size=int(np.floor(d0*pcent)))
        def D(j):
            return data["time"][j]-data["time"][j-1]
        Dfunc = np.vectorize(D)
        aux_arr = Dfunc(rand_i)
        return np.sum(aux_arr)/len(aux_arr)
        
    def rate5(self,xseed=20171224):
        data = self.lc
        np.random.seed(xseed)
        d0 = data.shape[0]
        pcent = 0.25
        rand_i = np.random.randint(1,high=d0-1,size=int(np.floor(d0*pcent)))
        F = lambda j: data["time"][j]-data["time"][j-1]
        aux_arr = map(F,rand_i)
        return np.sum(aux_arr)/len(aux_arr)
    
    def rate6(self,xseed=20171224):
        """the Fastest!"""
        data = self.lc
        np.random.seed(xseed)
        d0 = data.shape[0]
        pcent = 0.25
        rand_i = np.random.randint(1,high=d0-1,size=int(np.floor(d0*pcent)))
        F = lambda j: data["time"][j]-data["time"][j-1]
        Fvec = np.vectorize(F)
        aux_arr = Fvec(rand_i)
        return np.sum(aux_arr)/len(aux_arr)


if __name__ == "__main__":
    print "Starting"
    """Call the script using the path to time series
    Example: > python compare_loops.py lightcurves/src1/
    """
    folder = sys.argv[1]
    
    #           Fractal Dim and Lacunarity
    lacun = False
    DEPTH = 0
    for root,dirs,files in os.walk(folder):
        if root.count(os.sep) >= DEPTH:
            del dirs[:]
        for index,item in enumerate(files):
            fnm = os.path.join(root,item)
            if (".dat" in item) and (os.access(fnm,os.R_OK)):
                print fnm
                lc = LC.open_lc(fnm)
                TimeClass = Play(lc)
                tr1 = timeit.Timer(TimeClass.rate1).timeit(number=1000)
                tr2 = timeit.Timer(TimeClass.rate2).timeit(number=1000)
                tr3 = timeit.Timer(TimeClass.rate3).timeit(number=1000)
                tr4 = timeit.Timer(TimeClass.rate4).timeit(number=1000)
                tr5 = timeit.Timer(TimeClass.rate5).timeit(number=1000)
                tr6 = timeit.Timer(TimeClass.rate6).timeit(number=1000)
                print "Method rate1: {0} sec".format(tr1) 
                print "Method rate2: {0} sec".format(tr2) 
                print "Method rate3: {0} sec".format(tr3) 
                print "Method rate4: {0} sec".format(tr4) 
                print "Method rate5: {0} sec".format(tr5) 
                print "Method rate6: {0} sec".format(tr6) 
                
                print "exiting (debug)..."; exit()
            elif (".dat" in item) and ~(os.access(fnm,os.R_OK)):
                logging.error("No reading permission for {0}".format(fnm))

