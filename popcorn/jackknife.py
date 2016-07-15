from __future__ import division
from __future__ import print_function
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
              ' error using',self.nb,'blocks of size',self.bs,'.')
        self.psuedovalues, self.delete_values = self.jackknife(f,x,r)
        maskedps=np.ma.array(self.psuedovalues,mask=np.isnan(self.psuedovalues))
        nb = self.nb - np.isnan(self.psuedovalues).sum(0)
        self.cov = np.ma.cov(maskedps,rowvar=False)
        self.SE = np.sqrt(np.diag(self.cov)/nb)

    def get_blocks(self,N,args):
        bs0=200
        nb0 = 200
        minb = 10
        nblocks = np.min((np.max((minb,N/bs0)),nb0))
        bs = int(N/nblocks)
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
            sys.stdout.write("Jackknife iter: {0} \r".format(i+1))
            sys.stdout.flush()
            k = np.ones((x.shape[0]),dtype=np.bool)
            k[slice(*b)]=False
            xi = x[k]
            dv = np.asarray(f(xi))
            ps = self.nb*dv - (self.nb-1)*r
            psa.append(ps)
            dva.append(dv)
        sys.stdout.write("\n")
        return np.array(psa), np.array(dva)
