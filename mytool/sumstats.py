from __future__ import division
from __future__ import print_function
import numpy as np
import pandas as pd
import time
import sys
import argparse
import warnings
from IPython import embed

compliment = {'A':'T','T':'A','G':'C','C':'G',
              'a':'t','t':'a','g':'c','c':'g'}

class sumstats_1_trait(object):
    '''
    Class for Sumstats objects.
    '''
    def __init__(self,scores,args):
        data1, id_type = self.parse_input(args.sfile,args.old_format)
        data1 = data1.loc[~data1.index.duplicated()]
        print(len(data1),"SNPs detected in input.")
        self.id_type = id_type
        if self.id_type == 'pos':
            scores['id']='chr'+scores['chr'].map(str)+':'+scores['pos'].map(str)
        comm_snps = scores.index.intersection(data1.index)
        print(len(comm_snps),"SNPs remaining after merging with score file.")
        data1 = data1.loc[comm_snps]
        scores = scores.loc[comm_snps]
        align1=self.align_to_scores(data1,scores['a1'],scores['a2'],scores['af'])
        data = pd.DataFrame()
        #data = scores[['chr','pos','id','a1','a2','score']]
        data = scores.copy()
        data['N']=data1['N']
        data['Z']=align1*data1['Z']
        try:
            data['beta']=align1*data1['beta']
            data['SE'] = data1['SE']
        except KeyError:
            pass
        print(len(data),"SNPs remaining after filtering on AF and self-"
              "compliment SNPs.")
        self.data = data.loc[align1!=0]

    def parse_input(self,sfile,old_format):
        if not old_format:
            DF = pd.read_table(sfile,sep='\s*',engine='python')
            data = pd.DataFrame()
            try:
                data['id'] = DF['rsid']
                id_type = 'rsid'
            except KeyError:
                try:
                    data['id'] = 'chr'+DF['chr'].map(str)+':'+DF['pos'].map(str)
                    id_type = 'pos'
                except KeyError:
                    raise ValueError('Must provide either "rsid" or "chr" and "pos"')
            try:
                data['a1'] = DF['a1']
                data['a2'] = DF['a2']
            except KeyError:
                raise ValueError('Must provide allele names "a1" and "a2"')
            try:
                data['af'] = DF['af']
            except KeyError:
                data['af'] = np.nan
                warnings.warn('Warning: Study AF not provided. Defulating to'
                              'ignoring variants where minor allele is the '
                              'compliment of the major allele.')
            try:
                data['N'] = DF['N']
            except KeyError:
                raise ValueError('Must provide number of samples "N"')
            try:
                data['beta'] = DF['beta']
                data['SE'] = DF['SE']
                data['Z'] = data['beta']/data['SE']
            except KeyError:
                try:
                    data['Z'] = DF['Z']
                except KeyError:
                    raise ValueError('Must provide either signed Z-scores "Z" or'
                                     '"beta" and "SE"')
            data.index = data['id']
        else:
            data = pd.read_table(sfile,sep='\t',header=None,
                                 names=['chr','id','pos','af','a1','a2',
                                          'N1','N2','beta','SE','pv'])
            data['N'] = data['N1']+data['N2']
            data['Z'] = data['beta']/data['SE']
            data = data[['id','af','a1','a2','N','beta','SE','Z']]
            id_type='rsid'
        return data, id_type

    # Similar enough to the one in compute.py to be rolled into one function
    def align_to_scores(self, data, a1, a2, afr, tol=0.1):
        _data = data.copy()
        a1c = [compliment[a] for a in a1]
        a2c = [compliment[a] for a in a2]
        _data['a1c'] = [compliment[a] for a in _data['a1']]
        _data['a2c'] = [compliment[a] for a in _data['a2']]

        self_compliment=(_data['a2']==_data['a1c'])
        s = (_data['a1']==a1)&\
            (_data['a2']==a2)&\
            (~self_compliment)
        f = (_data['a1']==a2)&\
            (_data['a2']==a1)&\
            (~self_compliment)
        c = (_data['a1']==a1c)&\
            (_data['a2']==a2c)&\
            (~self_compliment)
        fac = (_data['a1']==a2c)&\
            (_data['a2']==a1c)&\
            (~self_compliment)
        af_s = abs(data['af']-afr) < tol
        af_f = abs(data['af']-(1-afr)) < tol
        saf = (_data['a1']==a1)&\
            (_data['a2']==a2)&\
            self_compliment & af_s
        faf = (_data['a1']==a2)&\
            (_data['a2']==a1)&\
            self_compliment & af_f
        caf = (_data['a1']==a1c)&\
            (_data['a2']==a2c)&\
            self_compliment & af_s
        facaf = (_data['a1']==a2c)&\
            (_data['a2']==a1c)&\
            self_compliment & af_f
        alignment = np.array(s+-1*f+c+-1*fac+saf+-1*faf+caf+-1*facaf)
        #keep = (alignment!=0)
        return alignment #data.loc[keep], alignment[keep]

class sumstats_2_trait(sumstats_1_trait):
    def __init__(self,scores,args):
        data1, id_type1 = self.parse_input(args.sfile1)
        data2, id_type2 = self.parse_input(args.sfile2)
        try:
            assert id_type1==id_type2
        except AssertionError:
            raise ValueError('Provided summary statistics files do not have'
                             ' the same SNP ID types.')
        print(len(data1),len(data2),"SNPs detected in input.")
        self.id_type = id_type1
        if self.id_type == 'pos':
            # Need to double check that this changes the index since id=index
            scores['id']='chr'+scores['chr'].map(str)+':'+scores['pos'].map(str)
        comm_snps = scores.index.intersection(data1.index).\
            intersection(data2.index)
        print(len(comm_snps),"SNPs remaining after merging with score file.")
        data1 = data1.loc[comm_snps]
        data2 = data2.loc[comm_snps]
        scores = scores.loc[comm_snps]
        try:
            af1, af2 = scores['af'], scores['af']
            self.two_pops = False
        except KeyError:
            af1, af2 = scores['af1'], scores['af2']
            self.two_pops = True
        align1=self.align_to_scores(data1,scores['a1'],scores['a2'],af1)
        align2=self.align_to_scores(data2,scores['a1'],scores['a2'],af2)
        # data = pd.DataFrame()
        # try:
        #     data = scores[['chr','pos','id','a1','a2','score']]
        # except KeyError:
        #     data = scores[['chr','pos','id','a1','a2',
        #                    'score1','score2','scoreX']]
        data = scores.copy()
        data['N1']=data1['N'] #Gives a warning I don't know why
        data['N2']=data2['N']
        data['Z1']=align1*data1['Z']
        data['Z2']=align2*data1['Z']
        try:
            data['beta1']=align1*data1['beta']
            data['beta2']=align2*data2['beta']
            data['SE1'] = data1['SE']
            data['SE2'] = data2['SE']
        except KeyError:
            pass
        try:
            # Haven't tested, probably doesn't work
            Ns = pd.read_table(args.Ns,index_col=0)
            comm = NS.index.intersection(data.index)
            data['Ns'] = Ns.loc[comm]
        except IOError:
            try:
                Ns = float(args.Ns)
                if Ns < 1.0:
                    data['Ns'] = Ns*data['N1']
                else:
                    data['Ns'] = Ns
            except TypeError:
                if self.two_pops == False:
                    warnings.warn('Sample overlap not specified,'
                                 ' assuming 0.')
                data['Ns'] = 0
        self.overlap = np.all(data['Ns'])
        if self.two_pops and self.overlap:
            raise ValueError('Sample overlap between two populations'
                             ' specified.')
        print(len(data),"SNPs remaining after filtering on AF and self-"
              "compliment SNPs.")

        self.data = data.loc[(align1!=0)&(align2!=0)]
