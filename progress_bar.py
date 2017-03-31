import os
import sys
import subprocess
import time
import numpy as np

class Toolbox():
    @classmethod
    def progress_bar(cls,iterator,Nposit,wait_time=0.25,bar_wide=30):
        '''Receives the actual iterator and the max number of items 
        Idea from: http://stackoverflow.com/questions/3002085/
        python-to-print-out-status-bar-and-percentage
        '''
        sys.stdout.write('\r')
        aux = (iterator*100/Nposit)-1
        sys.stdout.write('|{0:{1}}| {2}%'.format('='*iterator,bar_wide,aux))
        sys.stdout.flush()
        time.sleep(wait_time)
        return True
