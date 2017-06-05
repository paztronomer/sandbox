"""Simple script to copy a bunch of files living in different directories
"""
import os
import subprocess
import shlex
import argparse
import logging
import getpass
import numpy as np

class Toolbox():
    def __init__(self,keylock,input_lst,remote1,
                add_root=None,remote2=None,wildname=None):
        self.pw = keylock
        kw1 = {"dtype":"|S200","comments":"#"}
        if not (add_root is None):
            kw1["converters"] = {0: lambda x: os.path.join(add_root,x)}
        tmp = np.loadtxt(input_lst,**kw1)
        if tmp.size == 1:
            logging.error("Input list must have 2 or more entries")
            exit(1)
        if wildname is None:
            self.src = tmp
        else:
            self.src = list(map(lambda x:os.path.join(x,wildname),tmp))
        remote1 += ":"
        if not (remote2 is None):
            remote1 += remote2
        self.remo = remote1

    def one(self):
        """Do the copy
        """
        for i in self.src:
            print "\n\tCP-ing {0} to {1}".format(i,self.remo)
            cmd = "cp -v {0} {1}".format(i,self.remo)
            print cmd
            cmd = shlex.split(cmd)
            try:
                pB = subprocess.Popen(cmd,shell=False,
                                stdin=None,#subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True)
                outM,errM = pB.communicate()#input=self.pw+"\n")
                if errM:
                    print errM,"_"*30
                pB.wait()
            except:
                logging.error("Error in copying {0}".format(i))

if __name__=="__main__":
    h0 = "Code to copy a bunch of files living in different local"
    h0 += " directories to a single location"
    bun = argparse.ArgumentParser(description=h0)
    #positional
    h1 = "Filename (or complete path) to the list of source paths, from which"
    h1 += " data will be copied. Format of the list: one column with the path"
    h1 += " or filename"
    bun.add_argument("source_list",help=h1,type=str)
    h4 = "Directory in the remote location, where to save the files. If no"
    h4 += " value is given, files will be stored in the home"
    bun.add_argument("target",help=h4)
    #optional
    h3 = "Root path to be added as prefix, to the paths given in the input list"
    bun.add_argument("--root",help=h3,metavar="")
    h5 = "Wildcard for the files to be copied from source folders to target."
    h5 += " Example: *_hpix.fits"
    bun.add_argument("--wild",help=h5,metavar="")
    args = bun.parse_args()
    kw0 = vars(args)
    ls = [kw0["source_list"],kw0["target"]]
    di = dict()
    di["add_root"] = kw0["root"]
    di["wildname"] = kw0["wild"]
    #
    #x = getpass.getpass(prompt="Password to {0}:\n".format(kw0["host"]))
    #
    Toolbox(x,*ls,**di).one()
