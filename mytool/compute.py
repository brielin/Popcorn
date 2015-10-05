from __future__ import division
from __future__ import print_function
import numpy as np
import pandas as pd
import statsmodels.api as sm
import sys
from pysnptools.snpreader import Bed
from pysnptools.standardizer import Unit
from IPython import embed
from time import time

compliment = {'A':'T','T':'A','G':'C','C':'G',
              'a':'t','t':'a','g':'c','c':'g'}

class covariance_scores_1_pop(object):
    '''
    Class for storing covariance score objects and computing covariance scores

    Paramters
    ---------
    bfile: bed file name
    block_size: block size for computing scores
    block_type: type of block - SNP or BP

    Attributes
    ----------
    blocks: block boundaries
    block_type: type of block used to create the object
    N: sample size of panel used to determine scores
    M: number of SNPs
    scores: the covariance scores

    Methods
    -------
    get_blocks: determine the block boundaries
    compute: compute the covariance scores
    '''
    def __init__(self,args):
        if args.window_type not in ['BP','SNP']:
            raise ValueError('Window type not supported')
        bed_1 = Bed(args.bfile) #
        af1 = self.get_allele_frequency(bed_1,args) #
        snps_1 = (af1>args.maf)&(af1<1-args.maf) #
        if (args.from_bp is not None) and (args.to_bp is not None):
            k = (bed_1.pos[:,2]>args.from_bp)&(bed_1.pos[:,2]<args.to_bp)
            snps_1 = snps_1&k
        snps_to_use = bed_1.sid[snps_1]
        bed_1_index = np.sort(bed_1.sid_to_index(snps_to_use)) #
        pos = bed_1.pos[bed_1_index] #
        bim_1=pd.read_table(bed_1.filename+'.bim',header=None,
                            names=['chm','id','pos_mb','pos_bp','a1','a2'])
        af = af1[bed_1_index] #
        if args.afile is not None:
            a1 =  pd.read_table(args.afile,header=None,sep='\s*',
                                names=['id1','id2','theta'])
        else:
            a1 = None
        self.af = af
        self.M = len(bed_1_index) #
        self.windows = self.get_windows(pos,args) #
        self.chr = pos[:,0]
        self.pos = pos[:,2]
        self.id = bed_1.sid[bed_1_index]
        self.A1 = bim_1['a1'].loc[bed_1_index]
        self.A2 = bim_1['a2'].loc[bed_1_index]
        self.scores = self.compute(bed_1,bed_1_index,af,a1,args) #

    def get_windows(self,pos,args):
        if args.window_type == 'BP':
            coords = pos[:,2]
            ws = 1000*args.window_size
        elif args.window_type == 'SNP':
            coords = np.array(range(self.M))
            ws = args.window_size
        wl = []
        wr = []
        j=0
        for i in xrange(self.M):
            while j<self.M and abs(coords[j]-coords[i])>ws:
                j+=1
            wl.append(j)
        j=0
        for i in xrange(self.M):
            while j<self.M and wl[j] <= i:
                j += 1
            wr.append(j)
        return np.array([wl,wr]).T

    def _fast_var(self,X,m):
        return ((X - m)**2).mean(0)

    def get_allele_frequency(self,bed,args):
        s = args.SNPs_to_read
        af = np.zeros((bed.sid_count))
        var = np.zeros((bed.sid_count))
        if (args.from_bp is not None) and (args.to_bp is not None):
            k0 = np.where((bed.pos[:,2]>=args.from_bp))[0][0]
            k1 = np.where((bed.pos[:,2]<=args.to_bp))[0][-1]
            X = bed[:,k0:k1].read().val
            af[k0:k1] = X.mean(0)/2.0
            var[k0:k1] = self._fast_var(X,2*af[k0:k1])
        else:
            for i in xrange(int(np.ceil(bed.sid_count/s))):
                X = bed[:,i*s:(i+1)*s].read().val
                af[i*s:(i+1)*s] = X.mean(0)/2.0
                var[i*s:(i+1)*s] = self._fast_var(X,2*af[i*s:(i+1)*s])
        af[var==0]=0
        return af

    def _norm_data(self,X):
        m = np.nanmean(X,0)
        inds = np.where(np.isnan(X))
        X[inds]=np.take(m,inds[1])
        np.subtract(X,m,out=X)
        v = X.var(0)
        np.divide(X,np.sqrt(v),out=X)
        return X

    def _condition_ancestry(self,X,a):
        return sm.OLS(X,sm.add_constant(a)).fit().resid

    def compute(self,bed_1,bed_1_index,af,a1,args):
        N = bed_1.iid_count
        if args.per_allele:
            v1m = np.mean(2*af*(1-af))
            def func(a,i,j):
                af1 = af[i:j]
                v = (2*af1*(1-af1))/v1m
                v1j = 2*(af1[-1]*(1-af1[-1]))/v1m
                c = a**2
                return (v*c).sum(1), (v1j*(c)).sum(0)[0:-1]
        else:
            def func(a,i,j):
                c = a**2
                return c.sum(1), c.sum(0)[0:-1]
        t=time()
        scores = np.zeros((self.M))
        li,ri = self.windows[0]
        A1 = self._norm_data(bed_1[:,bed_1_index[li:ri]].read().val)
        if a1 is not None:
            A1 = self._condition_ancestry(A1,a1['theta'].values)
        R1 = np.dot(A1.T,A1/N)
        scores[li:ri] += func(R1,li,ri)[0]
        nstr = ri-li
        offset = 0
        #out1 = np.zeros((1,nstr-1))
        for i in xrange(ri,self.M,nstr):
            sys.stdout.write("SNP: %d, %f\r" % (i, time()-t))
            sys.stdout.flush()
            X1n= self._norm_data(bed_1[:,bed_1_index[i:(i+nstr)]].read().val)
            A1 = np.hstack((A1,X1n))
            if a1 is not None:
                A1 = self._condition_ancestry(A1,a1['theta'].values)
            for j in xrange(i,np.min((i+nstr,self.M))):
                lb,rb = self.windows[j]
                lbp = lb-offset
                jp = j-offset
                out1 = np.dot(np.atleast_2d(A1[:,jp]/N),A1[:,lbp:jp])
                # np.dot(np.atleast_2d(A1[:,jp]/N),A1[:,lbp:jp],out=out1)
                _out1 = np.hstack((out1,[[1]]))
                func_ret = func(_out1,lb,j+1)
                scores[lb:j] += func_ret[1]
                scores[j] += func_ret[0]
            if A1.shape[1] > args.SNPs_to_store:
                A1 = A1[:,nstr:]
                offset += nstr
        print(time()-t)
        return scores

    def write(self,args):
        f = open(args.out,'w')
        f.write('# M = '+str(self.M)+'\n')
        for l in zip(self.chr,self.pos,self.id,self.A1,self.A2,self.af,
                     self.scores):
            f.write('\t'.join(map(str,l))+'\n')

class covariance_scores_2_pop(covariance_scores_1_pop):
    def __init__(self,args):
        if args.window_type not in ['BP','SNP']:
            raise ValueError('Window type not supported')
        bed_1 = Bed(args.bfile1) #
        bed_2 = Bed(args.bfile2)
        af1 = self.get_allele_frequency(bed_1,args) #
        af2 = self.get_allele_frequency(bed_2,args)
        snps_1 = (af1>args.maf)&(af1<1-args.maf) #
        snps_2 = (af2>args.maf)&(af2<1-args.maf)
        if (args.from_bp is not None) and (args.to_bp is not None):
            k1 = (bed_1.pos[:,2]>args.from_bp)&(bed_1.pos[:,2]<args.to_bp)
            k2 = (bed_2.pos[:,2]>args.from_bp)&(bed_2.pos[:,2]<args.to_bp)
            snps_1 = snps_1&k1
            snps_2 = snps_2&k2
        snps_to_use = np.intersect1d(bed_1.sid[snps_1],bed_2.sid[snps_2])
        bed_1_index = np.sort(bed_1.sid_to_index(snps_to_use)) #
        bed_2_index = np.sort(bed_2.sid_to_index(snps_to_use))
        if not args.no_align:
            alignment,bed_1_index,bed_2_index =\
                self.align_alleles(bed_1,bed_1_index,af1,bed_2,bed_2_index,af2)
        else:
            alignment = np.ones(len(bed_1_index))
        pos = bed_1.pos[bed_1_index] #
        bim_1=pd.read_table(bed_1.filename+'.bim',header=None,
                            names=['chm','id','pos_mb','pos_bp','a1','a2'])
        af1 = af1[bed_1_index] #
        af2 = af2[bed_2_index]
        if args.afile1 is not None:
            a1 =  pd.read_table(args.afile,header=None,sep='\s*',
                                names=['id1','id2','theta'])
        else:
            a1 = None
        if args.afile2 is not None:
            a2 =  pd.read_table(args.afile,header=None,sep='\s*',
                                names=['id1','id2','theta'])
        else:
            a2 = None
        self.af1 = af1 #
        self.af2 = af2
        self.M = len(bed_1_index) #
        self.N = (bed_1.iid_count, bed_2.iid_count) #
        self.chr = pos[:,0]
        self.pos = pos[:,2]
        self.id = bed_1.sid[bed_1_index]
        self.A1 = bim_1['a1'].loc[bed_1_index]
        self.A2 = bim_1['a2'].loc[bed_1_index]
        self.windows = self.get_windows(pos,args) #
        self.scores1 = self.compute(bed_1,bed_1_index,af1,a1,args)
        self.scores2 = self.compute(bed_2,bed_2_index,af2,a2,args) #
        self.scoresX = self.compute2(bed_1,bed_1_index,bed_2,bed_2_index,
                                     alignment,a1,a2,args) #

    def write(self,args):
        f = open(args.out,'w')
        for l in zip(self.chr,self.pos,self.id,self.A1,self.A2,self.af1,
                     self.af2,self.scores1,self.scores2,self.scoresX):
            f.write('\t'.join(map(str,l))+'\n')

    # Align pop2 alleles relative to pop1 alleles
    def align_alleles(self,bed_1,bed_1_index,af1,bed_2,bed_2_index,af2,tol=0.1):
        bed_1_index = np.sort(bed_1_index)
        bed_2_index = np.sort(bed_2_index)
        af1 = af1[bed_1_index]
        af2 = af2[bed_2_index]
        bim_1=pd.read_table(bed_1.filename+'.bim',header=None,
                            names=['chm','id','pos_mb','pos_bp','a1','a2'])
        bim_1=bim_1.iloc[bed_1_index]
        bim_1.index = xrange(len(bed_1_index))
        bim_2=pd.read_table(bed_2.filename+'.bim',header=None,
                            names=['chm','id','pos_mb','pos_bp','a1','a2'])
        bim_2=bim_2.iloc[bed_2_index]
        bim_2.index = xrange(len(bed_2_index))
        bim_1['a1c'] = [compliment[a] for a in bim_1['a1']]
        bim_1['a2c'] = [compliment[a] for a in bim_1['a2']]
        bim_2['a1c'] = [compliment[a] for a in bim_2['a1']]
        bim_2['a2c'] = [compliment[a] for a in bim_2['a2']]

        self_compliment=(bim_1['a2']==bim_1['a1c'])
        s = (bim_1['a1']==bim_2['a1'])&\
            (bim_1['a2']==bim_2['a2'])&\
            (~self_compliment)
        f = (bim_1['a1']==bim_2['a2'])&\
            (bim_1['a2']==bim_2['a1'])&\
            (~self_compliment)
        c = (bim_1['a1']==bim_2['a1c'])&\
            (bim_1['a2']==bim_2['a2c'])&\
            (~self_compliment)
        fac = (bim_1['a1']==bim_2['a2c'])&\
            (bim_1['a2']==bim_2['a1c'])&\
            (~self_compliment)
        af_s = abs(af1-af2) < tol
        af_f = abs(af1-(1-af2)) < tol
        saf = (bim_1['a1']==bim_2['a1'])&\
            (bim_1['a2']==bim_2['a2'])&\
            self_compliment & af_s
        faf = (bim_1['a1']==bim_2['a2'])&\
            (bim_1['a2']==bim_2['a1'])&\
            self_compliment & af_f
        caf = (bim_1['a1']==bim_2['a1c'])&\
            (bim_1['a2']==bim_2['a2c'])&\
            self_compliment & af_s
        facaf = (bim_1['a1']==bim_2['a2c'])&\
            (bim_1['a2']==bim_2['a1c'])&\
            self_compliment & af_f
        alignment = np.array(s+-1*f+c+-1*fac+saf+-1*faf+caf+-1*facaf)
        keep = (alignment!=0)
        return alignment[keep], bed_1_index[keep], bed_2_index[keep]

    def compute2(self,bed_1,bed_1_index,bed_2,bed_2_index,alignment,a1,a2,args):
        if args.per_allele:
            v1m = np.mean(2*self.af1*(1-self.af1))
            v2m = np.mean(2*self.af2*(1-self.af2))
            def func(a,b,i,j):
                af1 = self.af1[i:j]
                af2 = self.af2[i:j]
                v = np.sqrt(2*af1*(1-af1)*2*af2*(1-af2))/np.sqrt(v1m*v2m)
                v1j = 2*(af1[-1]*(1-af1[-1]))/v1m
                v2j = 2*(af2[-1]*(1-af2[-1]))/v2m
                c = a*b
                return (v*c).sum(1), (np.sqrt(v1j*v2j)*(c)).sum(0)[0:-1]
        else:
            def func(a,b,i,j):
                c = a*b
                return c.sum(1), c.sum(0)[0:-1]
        t=time()
        scores = np.zeros((self.M))
        li,ri = self.windows[0]
        A1 = bed_1[:,bed_1_index[li:ri]].read().standardize(Unit()).val
        A2 = bed_2[:,bed_2_index[li:ri]].read().standardize(Unit()).val
        if a1 is not None:
            A1 = self._condition_ancestry(A1,a1['theta'].values)
        if a2 is not None:
            A2 = self._condition_ancestry(A2,a2['theta'].values)
        R1 = np.dot(A1.T,A1/self.N[0])
        R2 = np.dot(A2.T,A2/self.N[1])
        align_mat = np.outer(alignment[li:ri],alignment[li:ri])
        scores[li:ri] += func(R1,align_mat*R2,li,ri)[0]
        nstr = ri-li
        offset = 0
        #out1 = np.zeros((1,nstr-1))
        #out2 = np.zeros((1,nstr-1))
        for i in xrange(ri,self.M,nstr):
            sys.stdout.write("SNP: %d, %f\r" % (i,time()-t))
            sys.stdout.flush()
            X1n= bed_1[:,bed_1_index[i:(i+nstr)]].read().standardize(Unit()).val
            X2n= bed_2[:,bed_2_index[i:(i+nstr)]].read().standardize(Unit()).val
            A1 = np.hstack((A1,X1n))
            A2 = np.hstack((A2,X2n))
            if a1 is not None:
                A1 = self._condition_ancestry(A1,a1['theta'].values)
            if a2 is not None:
                A2 = self._condition_ancestry(A2,a2['theta'].values)
            for j in xrange(i,np.min((i+nstr,self.M))):
                lb,rb = self.windows[j]
                lbp = lb-offset
                jp = j-offset
                # np.dot(np.atleast_2d(A1[:,jp]/self.N[0]),A1[:,lbp:jp],out=out1)
                # np.dot(np.atleast_2d(A2[:,jp]/self.N[1]),A2[:,lbp:jp],out=out2)
                out1=np.dot(np.atleast_2d(A1[:,jp]/self.N[0]),A1[:,lbp:jp])
                out2=np.dot(np.atleast_2d(A2[:,jp]/self.N[1]),A2[:,lbp:jp])
                align_mat = np.atleast_2d(alignment[j]*alignment[lb:j])
                _out1 = np.hstack((out1,[[1]]))
                _out2 = np.hstack((align_mat*out2,[[1]]))
                func_ret = func(_out1,_out2,lb,j+1)
                scores[lb:j] += func_ret[1]
                scores[j] += func_ret[0]
            if A1.shape[1] > args.SNPs_to_store:
                A1 = A1[:,nstr:]
                A2 = A2[:,nstr:]
                offset += nstr
        print(time()-t)
        return scores
