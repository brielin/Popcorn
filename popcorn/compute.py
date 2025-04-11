import numpy as np
import pandas as pd
import bottleneck as bn
import sys
from pysnptools.snpreader import Bed
from pysnptools.standardizer import Unit
# from IPython import embed
from time import time

compliment = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'g':'c', 'c':'g', 1:3, 2:4}

class covariance_scores_1_pop(object):
    '''
    Class for storing covariance score objects and computing covariance scores
    Paramters
    ---------
    bfile: bed file name
    block_size: block size for computing scores
    block_type: type of block - SNP or KBP

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
        if args.window_type not in ['KBP','SNP']:
            raise ValueError('Window type not supported')
        bed_1 = Bed(args.bfile,count_A1=False) #
        af1 = self.get_allele_frequency(bed_1,args) #
        print(len(af1), "SNPs in file 1")
        snps_1 = (af1>args.maf)&(af1<1-args.maf) #
        print(np.sum(snps_1), "SNPs in file 1 after MAF filter")
        # Omit SNPs with NA values for h2weight
        if args.h2weight:
            data = pd.read_table(args.bfile+'.h2weight',header=None,
                                 names=['SNP','h2weight'],
                                 index_col=False)
            if(len(data['SNP']) != len(bed_1.sid) or
                       (data['SNP'] == bed_1.sid).min() == False):
                raise ValueError('SNPs disagree between bed/bim/fam and h2weight files')
            h2weight = data['h2weight']
            snps_1 = snps_1 & ~h2weight.isnull().values
            print(np.sum(snps_1), "SNPs in file 1 after extracting non-NA h2weight")
            del(data)
        if (args.from_bp is not None) and (args.to_bp is not None):
            k = (bed_1.pos[:,2]>args.from_bp)&(bed_1.pos[:,2]<args.to_bp)
            snps_1 = snps_1&k
        snps_to_use = bed_1.sid[snps_1]
        if args.extract is not None:
            keep = np.array([l.strip() for l in open(args.extract,'r')])
            snps_to_use = np.intersect1d(snps_to_use,keep)
            print(len(snps_to_use),"SNPs remaining after extraction")
        bed_1_index = np.sort(bed_1.sid_to_index(snps_to_use)) #
        pos = bed_1.pos[bed_1_index] #
        bim1_filename = bed_1.filename
        bim_1 = pd.read_table(bim1_filename.replace('.bed','')+'.bim',header=None,
                            names=['chm','id','pos_mb','pos_bp','a1','a2'])
        af = af1[bed_1_index] #
        # if args.afile is not None:
        #     a1 =  pd.read_table(args.afile,header=None,sep='\s*',
        #                         names=['id1','id2','theta'])
        # else:
        a1 = None
        try:
            h2weight = h2weight[bed_1_index].values
        except NameError:
            h2weight = None
        self.af = af
        self.M = len(bed_1_index) #
        self.windows = self.get_windows(pos,args) #
        self.chr = pos[:,0]
        self.pos = pos[:,2]
        self.id = bed_1.sid[bed_1_index]
        self.A1 = bim_1['a1'].loc[bed_1_index]
        self.A2 = bim_1['a2'].loc[bed_1_index]
        self.scores = self.compute(bed_1,bed_1_index,af,a1,h2weight,args) #

    def get_windows(self,pos,args):
        if args.window_type == 'KBP':
            coords = pos[:,2]
            ws = 1000*args.window_size
        elif args.window_type == 'SNP':
            coords = np.array(range(self.M))
            ws = args.window_size
        wl = []
        wr = []
        j=0
        for i in range(self.M):
            while j<self.M and abs(coords[j]-coords[i])>ws:
                j+=1
            wl.append(j)
        j=0
        for i in range(self.M):
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
            af[k0:k1] = bn.nanmean(X,0)/2.0
            var[k0:k1] = self._fast_var(X,2*af[k0:k1])
        else:
            for i in range(int(np.ceil(bed.sid_count/s))):
                X = bed[:,i*s:(i+1)*s].read().val
                af[i*s:(i+1)*s] = bn.nanmean(X,0)/2.0
                var[i*s:(i+1)*s] = self._fast_var(X,2*af[i*s:(i+1)*s])
        af[var==0]=0
        return af

    # replace missing genotypes with mean
    def _norm_data(self,X):
        m = bn.nanmean(X,0)
        inds = np.where(np.isnan(X))
        X[inds]=np.take(m,inds[1])
        np.subtract(X,m,out=X)
        v = X.var(0)
        np.divide(X,np.sqrt(v),out=X)
        return X

    # TODO: test ancestry conditioning
    # def _condition_ancestry(self,X,a):
    #     return sm.OLS(X,sm.add_constant(a)).fit().resid

    def compute(self,bed_1,bed_1_index,af,a1,h2weight,args):
        N = bed_1.iid_count
        if args.gen_effect:
            v1m = np.mean(2*af*(1-af))
            def func(a, i, j):
                af1 = af[i:j]
                v = (2*af1*(1-af1))/v1m
                v1j = 2*(af1[-1]*(1-af1[-1]))/v1m
                c = a**2
                if not args.use_bias:
                    c = c - (1-c)/(N-2)
                return (v*c).sum(1), (v1j*(c)).sum(0)[0:-1]
        elif args.h2weight:
            v1m = np.mean(h2weight)
            def func(a, i, j):
                h2w1 = h2weight[i:j]
                v = h2w1/v1m
                v1j = h2w1[-1]/v1m
                c = a**2
                if not args.use_bias:
                    c = c - (1-c)/(N-2)
                return (v*c).sum(1), (v1j*(c)).sum(0)[0:-1]
        else:
            def func(a, i, j):
                c = a**2
                if not args.use_bias:
                    c = c - (1-c)/(N-2)
                return c.sum(1), c.sum(0)[0:-1]
        t=time()
        scores = np.zeros((self.M))
        li,ri = self.windows[0]
        A1 = self._norm_data(bed_1[:,bed_1_index[li:ri]].read().val)
        # if a1 is not None:
        #    A1 = self._condition_ancestry(A1,a1['theta'].values)
        R1 = np.dot(A1.T,A1/N)
        scores[li:ri] += func(R1,li,ri)[0]
        nstr = np.max((1000,ri-li))
        #nstr = ri-li
        offset = 0
        #out1 = np.zeros((1,nstr-1))
        for i in range(ri,self.M,nstr):
            sys.stdout.write("SNP: %d, %f\r" % (i, time()-t))
            sys.stdout.flush()
            X1n= self._norm_data(bed_1[:,bed_1_index[i:(i+nstr)]].read().val)
            A1 = np.hstack((A1,X1n))
            # if a1 is not None:
            #     A1 = self._condition_ancestry(A1,a1['theta'].values)
            for j in range(i,np.min((i+nstr,self.M))):
                lb,rb = self.windows[j]
                lbp = lb-offset
                jp = j-offset
                out1 = np.dot(np.atleast_2d(A1[:,jp]/N),A1[:,lbp:jp])
                # np.dot(np.atleast_2d(A1[:,jp]/N),A1[:,lbp:jp],out=out1)
                _out1 = np.hstack((out1,[[1]]))
                func_ret = func(_out1,lb,j+1)
                try:
                    scores[lb:j] += func_ret[1]
                except ValueError:
                    raise ValueError("Error when setting scores."\
                          " Block width may exceed number of SNPs being stored"\
                          " Try increasing --SNPs_to_store"\)
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
        if args.window_type not in ['KBP','SNP']:
            raise ValueError('Window type not supported')
        # Open files
        bed_1 = Bed(args.bfile1,count_A1=False) #
        bed_2 = Bed(args.bfile2,count_A1=False)
        # Get indel locations, if any
        bim1_filename = bed_1.filename
        bim2_filename = bed_2.filename
        bim_1=pd.read_table(bim1_filename.replace('.bed','')+'.bim',header=None,
                            names=['chm','id','pos_mb','pos_bp','a1','a2'])
        bim_2=pd.read_table(bim2_filename.replace('.bed','')+'.bim',header=None,
                            names=['chm','id','pos_mb','pos_bp','a1','a2'])
        is_indel_1 = np.array([(len(str(a1))>1)|(len(str(a2))>1)  for a1,a2 in bim_1[['a1','a2']].values])
        is_indel_2 = np.array([(len(str(a1))>1)|(len(str(a2))>1)  for a1,a2 in bim_2[['a1','a2']].values])
        # Make sure two SNPs don't have the same position
        is_duplicated_bp_1=bim_1.pos_bp.duplicated()
        is_duplicated_bp_2=bim_2.pos_bp.duplicated()
        # Make sure two SNPs don't have the same ID
        is_duplicated_id_1=bim_1.id.duplicated()
        is_duplicated_id_2=bim_2.id.duplicated()
        # Get allele frequencies
        af1 = self.get_allele_frequency(bed_1,args) #
        af2 = self.get_allele_frequency(bed_2,args)

        print(len(af1), "Variants in file 1")
        print(len(af2), "Variants in file 2")

        # Get good SNPs
        snps_1 = (af1>args.maf)&(af1<1-args.maf)&(~is_indel_1)&(~is_duplicated_bp_1)&(~is_duplicated_id_1) #
        snps_2 = (af2>args.maf)&(af2<1-args.maf)&(~is_indel_2)&(~is_duplicated_bp_2)&(~is_duplicated_id_2)
        print(np.sum(snps_1), "SNPs in file 1 after MAF and indel filter")
        print(np.sum(snps_2), "SNPs in file 2 after MAF and indel filter")
        # Omit SNPs with NA values for h2weight
        if (args.h2weight):
            data = pd.read_table(args.bfile1+'.h2weight',header=None,
                                 names=['SNP','h2weight'],
                                 index_col=False)
            if(len(data['SNP']) != len(bed_1.sid) or
                       (data['SNP'] == bed_1.sid).min() == False):
                raise ValueError('SNPs disagree between bed/bim/fam and h2weight files 1')
            h2weight_1 = data['h2weight']
            snps_1 = snps_1 & ~h2weight_1.isnull().values
            print(np.sum(snps_1), "SNPs in file 1 after extracting non-NA h2weight")
            del(data)
            data = pd.read_table(args.bfile2+'.h2weight',header=None,
                                 names=['SNP','h2weight'],
                                 index_col=False)
            if(len(data['SNP']) != len(bed_2.sid) or
                       (data['SNP'] == bed_2.sid).min() == False):
                raise ValueError('SNPs disagree between bed/bim/fam and h2weight files 2')
            h2weight_2 = data['h2weight']
            snps_2 = snps_2 & ~h2weight_2.isnull().values
            print(np.sum(snps_2), "SNPs in file 2 after extracting non-NA h2weight")
            del(data)
        if (args.from_bp is not None) and (args.to_bp is not None):
            k1 = (bed_1.pos[:,2]>args.from_bp)&(bed_1.pos[:,2]<args.to_bp)
            k2 = (bed_2.pos[:,2]>args.from_bp)&(bed_2.pos[:,2]<args.to_bp)
            snps_1 = snps_1&k1
            snps_2 = snps_2&k2
        snps_to_use = np.intersect1d(bed_1.sid[snps_1.values],bed_2.sid[snps_2.values])
        print(len(snps_to_use),"SNPs common in both populations")
        if args.extract is not None:
            keep = np.array([l.strip() for l in open(args.extract,'r')])
            print(len(keep),"SNPs to extract")
            snps_to_use = np.intersect1d(snps_to_use,keep)
            print(len(snps_to_use),"SNPs remaining after extraction")
        bed_1_index = np.sort(bed_1.sid_to_index(snps_to_use)) #
        bed_2_index = np.sort(bed_2.sid_to_index(snps_to_use))
        if not args.no_align:
            alignment,bed_1_index,bed_2_index =\
                self.align_alleles(bed_1,bed_1_index,af1,bed_2,bed_2_index,af2)
        else:
            alignment = np.ones(len(bed_1_index))

        pos = bed_1.pos[bed_1_index] #
        af1 = af1[bed_1_index] #
        af2 = af2[bed_2_index]
        af2[alignment==-1] = (1-af2[alignment==-1])
        # if args.afile1 is not None:
        #     a1 =  pd.read_table(args.afile,header=None,sep='\s*',
        #                         names=['id1','id2','theta'])
        # else:
        a1 = None
        # if args.afile2 is not None:
        #     a2 =  pd.read_table(args.afile,header=None,sep='\s*',
        #                         names=['id1','id2','theta'])
        # else:
        a2 = None
        try:
            h2weight_1 = h2weight_1[bed_1_index].values
            h2weight_2 = h2weight_2[bed_2_index].values
        except NameError:
            h2weight_1 = None
            h2weight_2 = None
        self.af1 = af1 #
        self.af2 = af2
        self.h2weight_1 = h2weight_1
        self.h2weight_2 = h2weight_2
        self.M = len(bed_1_index) #
        self.N = (bed_1.iid_count, bed_2.iid_count) #
        self.chr = pos[:,0]
        self.pos = pos[:,2]
        self.id = bed_1.sid[bed_1_index]
        self.A1 = bim_1['a1'].iloc[bed_1_index]
        self.A2 = bim_1['a2'].iloc[bed_1_index]
        self.windows = self.get_windows(pos,args) #
        self.scores1 = self.compute(bed_1,bed_1_index,af1,a1,h2weight_1,args)
        self.scores2 = self.compute(bed_2,bed_2_index,af2,a2,h2weight_2,args) #
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
        bim1_filename = bed_1.filename
        bim_1=pd.read_table(bim1_filename.replace('.bed','')+'.bim',header=None,
                            names=['chm','id','pos_mb','pos_bp','a1','a2'])
        bim_1=bim_1.iloc[bed_1_index]
        bim_1.index = range(len(bed_1_index))
        bim2_filename = bed_2.filename
        bim_2=pd.read_table(bim2_filename.replace('.bed','')+'.bim',header=None,
                            names=['chm','id','pos_mb','pos_bp','a1','a2'])
        bim_2=bim_2.iloc[bed_2_index]
        bim_2.index = range(len(bed_2_index))
        a11c = []
        a12c = []
        a21c = []
        a22c = []
        for a11,a12,a21,a22 in \
                pd.concat([bim_1[['a1','a2']],bim_2[['a1','a2']]],axis=1).values:
            try:
                a11c.append(compliment[a11])
                a12c.append(compliment[a12])
                a21c.append(compliment[a21])
                a22c.append(compliment[a22])
            except KeyError:
                a11c.append('N')
                a12c.append('N')
                a21c.append('N')
                a22c.append('N')
        non_snp = (np.array(a11c,dtype='str') == 'N')
        bim_1['a1c'] = a11c
        bim_1['a2c'] = a12c
        bim_2['a1c'] = a21c
        bim_2['a2c'] = a22c
        # bim_1['a1c'] = [compliment[a] for a in bim_1['a1']]
        # bim_1['a2c'] = [compliment[a] for a in bim_1['a2']]
        # bim_2['a1c'] = [compliment[a] for a in bim_2['a1']]
        # bim_2['a2c'] = [compliment[a] for a in bim_2['a2']]
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
        keep = (alignment!=0)&(~non_snp)
        return alignment[keep], bed_1_index[keep], bed_2_index[keep]

    def compute2(self,bed_1,bed_1_index,bed_2,bed_2_index,alignment,a1,a2,args):
        if args.gen_effect:
            v1m = np.mean(2*self.af1*(1-self.af1))
            v2m = np.mean(2*self.af2*(1-self.af2))
            def func(a,b,i,j):
                af1 = self.af1[i:j]
                af2 = self.af2[i:j]
                v = np.sqrt(2*af1*(1-af1)*2*af2*(1-af2))/np.sqrt(v1m*v2m)
                v1j = 2*(af1[-1]*(1-af1[-1]))/v1m
                v2j = 2*(af2[-1]*(1-af2[-1]))/v2m
                c = a*b
                if not args.use_bias:
                    correction = (1 - np.abs(c))/np.sqrt((self.N[0]-2)*(self.N[1]-2))
                    c = c - np.sign(c)*correction
                return (v*c).sum(1), (np.sqrt(v1j*v2j)*(c)).sum(0)[0:-1]
        elif args.h2weight:
            v1m = np.mean(self.h2weight_1)
            v2m = np.mean(self.h2weight_2)
            def func(a,b,i,j):
                h2w1 = self.h2weight_1[i:j]
                h2w2 = self.h2weight_2[i:j]
                v = np.sqrt(h2w1*h2w2)/np.sqrt(v1m*v2m)
                v1j = h2w1[-1]/v1m
                v2j = h2w2[-1]/v2m
                c = a*b
                if not args.use_bias:
                    correction = (1 - np.abs(c))/np.sqrt((self.N[0]-2)*(self.N[1]-2))
                    c = c - np.sign(c)*correction
                return (v*c).sum(1), (np.sqrt(v1j*v2j)*(c)).sum(0)[0:-1]
        else:
            def func(a,b,i,j):
                c = a*b
                if not args.use_bias:
                    correction = (1 - np.abs(c))/np.sqrt((self.N[0]-2)*(self.N[1]-2))
                    c = c - np.sign(c)*correction
                return c.sum(1), c.sum(0)[0:-1]
        t=time()
        scores = np.zeros((self.M))
        li,ri = self.windows[0]
        A1 = bed_1[:,bed_1_index[li:ri]].read().standardize(Unit()).val
        A2 = bed_2[:,bed_2_index[li:ri]].read().standardize(Unit()).val
        # if a1 is not None:
        #     A1 = self._condition_ancestry(A1,a1['theta'].values)
        # if a2 is not None:
        #     A2 = self._condition_ancestry(A2,a2['theta'].values)
        R1 = np.dot(A1.T,A1/self.N[0])
        R2 = np.dot(A2.T,A2/self.N[1])
        align_mat = np.outer(alignment[li:ri],alignment[li:ri])
        scores[li:ri] += func(R1,align_mat*R2,li,ri)[0]
        nstr = np.max((1000,ri-li))
        #nstr = ri-li
        offset = 0
        #out1 = np.zeros((1,nstr-1))
        #out2 = np.zeros((1,nstr-1))
        for i in range(ri,self.M,nstr):
            sys.stdout.write("SNP: %d, %f\r" % (i,time()-t))
            sys.stdout.flush()
            X1n= bed_1[:,bed_1_index[i:(i+nstr)]].read().standardize(Unit()).val
            X2n= bed_2[:,bed_2_index[i:(i+nstr)]].read().standardize(Unit()).val
            A1 = np.hstack((A1,X1n))
            A2 = np.hstack((A2,X2n))
            # if a1 is not None:
            #     A1 = self._condition_ancestry(A1,a1['theta'].values)
            # if a2 is not None:
            #     A2 = self._condition_ancestry(A2,a2['theta'].values)
            for j in range(i,np.min((i+nstr,self.M))):
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
