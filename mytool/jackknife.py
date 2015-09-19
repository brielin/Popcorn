from __future__ import division
from __future__ import print_function
#from my_tool import sumstats
import numpy as np
import pandas as pd
import time
import sys
import argparse
import logging
import traceback
from scipy import optimize
from IPython import embed
from time import time

class jackknife(object):
    def __init__(self,f,x,r,args):
        N = len(x)
        self.blocks, self.bs = self.get_blocks(N,args)
        self.nb = len(self.blocks)
        print('Performing jackknife estimation of the standard'
              'error using',self.nb,'blocks of size',self.bs,'.')
        self.psuedovalues, self.delete_values = self.jackknife(f,x,r)
        self.cov = np.cov(self.psuedovalues,rowvar=False)
        self.SE = np.sqrt(np.diag(self.cov)/self.nb)

    def get_blocks(self,N,args):
        if N < 1000:
            bs = N/10
            nblocks = N/bs
        elif N < 10000:
            bs = 100
            nblocks = N/bs
        else:
            nblocks=200
            bs = N/nblocks
        A = np.floor(np.linspace(0,N,nblocks+1)).astype(int)
        return zip(A[:-1],A[1:]), bs

    def jackknife(self,f,x,r):
        psa = []
        dva = []
        r = np.asarray(r)
        for i,b in enumerate(self.blocks):
            # This is not the fastest method for np arrays
            # but it appears to be the fastest if x is a pandas
            # dataframe, plus the syntax is the same
            sys.stdout.write("Jackknife iter: %d \r" % i)
            sys.stdout.flush()
            k = np.ones((x.shape[0]),dtype=np.bool)
            k[slice(*b)]=False
            xi = x[k]
            dv = np.asarray(f(xi))
            ps = self.nb*dv - (self.nb-1)*r
            psa.append(ps)
            dva.append(dv)
        return np.array(psa), np.array(dva)
