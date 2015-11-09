from __future__ import division
from __future__ import print_function
import numpy as np
import pandas as pd
import sys
import argparse
import logging
import traceback
from scipy import optimize, stats
from IPython import embed
from time import time
from . import jackknife

class fit_by_region(object):
    def __init__(self,data,scores,args,t='h1'):
        types = {'h1':fit_h1,'pg':fit_pg,'pg_pe':fit_pg_pe}
        regions = pd.read_table(args.regions,sep='\s*')
        if regions['chr'][0][0]=='c':
            regions['chr'] = map(lambda x: int(x[3:]), regions['chr'])
        regions['M'] = self.get_M_by_region(regions,scores)
        fit_list = []
        for chm,start,end,M in regions.values:
            keep = (data['chr']==chm)&(data['pos']>start)&(data['pos']<end)
            datak = data[keep]
            fitk = types[t](datak,args,M=M)
            fit_list.append(fitk)
        region_index = regions['chr'].map(str)+':'+regions['start'].map(str)+\
            ':'+regions['stop'].map(str)
        self.fit_list = fit_list
        try:
            self.res_h = pd.concat(map(lambda x: x.res.loc['h'], fit_list),
                                   axis=1, ignore_index=True).T
            self.res_h.index = region_index
        except KeyError:
            pass
        try:
            self.res_h1 = pd.concat(map(lambda x: x.res.loc['h1'], fit_list),
                                    axis=1, ignore_index=True).T
            self.res_h1.index = region_index
            self.res_h2 = pd.concat(map(lambda x: x.res.loc['h2'], fit_list),
                                    axis=1, ignore_index=True).T
            self.res_h2.index = region_index
            self.res_pg = pd.concat(map(lambda x: x.res.loc['pg'], fit_list),
                                    axis=1, ignore_index=True).T
            self.res_pg.index = region_index
            self.res_pe = pd.concat(map(lambda x: x.res.loc['pe'], fit_list),
                                    axis=1, ignore_index=True).T
            self.res_pe.index = region_index
        except:
            pass

    def get_M_by_region(self,regions,scores):
        M = []
        for chm,start,end in regions.values:
            Mi = ((scores['chr']==chm)&(scores['pos']>start)&(scores['pos']<end)).sum()
            M.append(Mi)
        return np.array(M)

    def write(self,outfile):
        for k in self.__dict__.keys():
            try:
                self.__dict__[k].to_csv(outfile+'.'+k,sep='\t')
                print(self.__dict__[k].to_string())
            except AttributeError:
                pass


class fit_h1(object):
    def __init__(self,data,args,jk=True,M=None):
        if M is None: self.M = data.shape[0]
        else: self.M = M
        self.h_res, self.sy, self.ll = self.__call__(data,args)
        #print(self.h_res.x)
        self.null_ll = np.array([self.ll([1.0,0.0]), 0.0])
        #self.null_ll = np.array([self.ll(0.0), 0.0])
        self.alt_ll = np.array([self.h_res.fun, 0.0])
        res = np.array([self.h_res.x[1], self.sy])
        #res = np.array([self.h_res.x, self.sy])
        self.res = pd.DataFrame(np.vstack((res,np.tile(np.nan,len(res)))).T,
                                index=['h','sy'],columns=['Val','SE'])
        def close_call(x):
            res = self.__call__(x,args)
            return res[0].x[1], res[1]
            #return res[0].x, res[1]
        if jk and not args.no_jackknife:
            self.jackknife = jackknife.jackknife(close_call, data, res, args)
            self.res['SE'] = self.jackknife.SE
        else:
            self.jackknife = None
        if (args.K1 is not None)&(args.P1 is not None):
            self.convert_to_liability(args.K1,args.P1)
        self.res['LR'] = 2*(self.null_ll-self.alt_ll)
        self.res['Z'] = self.res['Val']/self.res['SE']
        self.res['P (LRT)'] = 1-stats.chi2.cdf(self.res['LR'],1)
        self.res['P (Z)'] = 1-stats.chi2.cdf(self.res['Z']**2,1)

    def __call__(self,data,args):
        W = 1.0/data['score']
        f=lambda x: self.nll(x,data['Z'],data['N'],self.M,data['score'],W=W)
        f2=lambda x: self.nll_no_intercept(x,data['Z'],data['N'],self.M,data['score'])
        # h = optimize.minimize_scalar(f,bounds=(0.0,1.0),method='bounded',
        #                              tol=args.tol, options={'disp':args.v})
        # sy = self.estimate_sy(data,h.x,data.shape[0])
        for x1 in np.arange(0.01,1,0.05):
            x0=[1.0,x1]
            h = optimize.minimize(f,x0,bounds=((None,None),(0.0,1.0)),
                                  options={'disp':args.v})
            if not h.success:
                continue
            else:
                break
        if not h.success:
            sys.stderr.write(h.message+'\n')
            raise ValueError
        sy = self.estimate_sy(data,h.x[1],data.shape[0])
        return h, sy, f

    def estimate_sy(self,data,h1,M):
        try:
            B = data['beta']
        except KeyError:
            return np.nan
        var = 2*data['af']*(1-data['af'])
        A = (1/(data['N']*var))*(1 + (data['N']/M)*h1*data['score'])
        return np.sum(B**2/A)/M

    def nll(self,x,Z,N,M,score,W=None):
        c,h = x
        S = (c + (N/M)*h*score)
#        S = (1 + (N/M)*x*score)
        if W is None:
            l = -0.5*M*np.log(2*np.pi) - 0.5*np.sum(np.log(S))\
                - 0.5*np.sum(Z**2/S)
        else:
            l = -0.5*M*np.log(2*np.pi) - 0.5*np.sum(W*np.log(S))\
                - 0.5*np.sum(W*(Z**2/S))
        return -1.0*l

    def nll_no_intercept(self,x,Z,N,M,score):
        S = (1 + (N/M)*x*score)
        l = -0.5*M*np.log(2*np.pi) - 0.5*np.sum(np.log(S))\
            - 0.5*np.sum(Z**2/S)
        return -1.0*l

    def write(self,outfile):
        print(self.res.to_string())
        self.res.to_csv(outfile,sep='\t',na_rep='NaN')

    def convert_to_liability(self,K,P):
        self.res['Lia']=np.nan
        ho=self.res['Val']['h']
        tau=stats.norm.ppf(1-K)
        hl=ho*((K*(1-K))**2)/(stats.norm.pdf(tau)**2*P*(1-P))
        self.res['Lia']['h']=hl
        return hl

class fit_pg(fit_h1):
    def __init__(self,data,args,M=None):
        if M is None: self.M = data.shape[0]
        else: self.M = M
        self.h1_res, self.sy1, self.h2_res, self.sy2, self.pg_res, self.ll =\
            self.__call__(data,args)
        h1r = self.h1_res.h_res.x[1]
        h2r = self.h2_res.h_res.x[1]
        pgr = self.pg_res.x
        hgr = pgr*np.sqrt(h1r*h2r)
        res = np.array([h1r, self.sy1, h2r, self.sy2, hgr, pgr])
        self.null_ll = np.array([self.h1_res.null_ll[0], 0.0,
                                 self.h2_res.null_ll[0], 0.0,self.ll(0.0),self.ll(0.0)])
        self.alt_ll = np.array([self.h1_res.alt_ll[0], 0.0,
                                self.h2_res.alt_ll[0], 0.0,self.pg_res.fun,self.pg_res.fun])
        self.res = pd.DataFrame(np.vstack((res,np.tile(np.nan,len(res)))).T,
                                index=['h1','sy1','h2','sy2','hg','pg'],
                                columns=['Val','SE'])
        print(self.res)
        def close_call(x):
            res = self.__call__(x,args)
            hg = res[4].x*np.sqrt(res[0].h_res.x[1]*res[2].h_res.x[1])
            return np.array([res[0].h_res.x[1], res[1], res[2].h_res.x[1],
                             res[3], hg, res[4].x])
        if not args.no_jackknife:
            self.jackknife = jackknife.jackknife(close_call,data,res,args)
            self.res['SE'] = self.jackknife.SE
        else:
            self.jackknife = None
        if (args.K1 is not None)&(args.P1 is not None)\
                &(args.K2 is not None)&(args.P2 is not None):
            self.convert_to_liability(args.K1,args.P1,args.K2,args.P2)
        self.res['LR'] = 2*(self.null_ll-self.alt_ll)
        self.res['Z'] = self.res['Val']/self.res['SE']
        self.res['P (LRT)'] = 1-stats.chi2.cdf(self.res['LR'],1)
        self.res['P (Z)'] = 1-stats.chi2.cdf(self.res['Z']**2,1)

    def __call__(self,data,args):
        data1 = data.copy()
        data2 = data.copy()
        try: #two populations
            data1[['N','af','Z','score']] = data[['N1','af1','Z1','score1']]
            data2[['N','af','Z','score']] = data[['N2','af2','Z2','score2']]
        except KeyError: #one population
            data1[['N','Z']] = data[['N1','Z1']]
            data2[['N','Z']] = data[['N2','Z2']]
        h1 = fit_h1(data1,args,jk=False,M=self.M)
        h2 = fit_h1(data2,args,jk=False,M=self.M)
        try:
            _fts = data['score1'] # A bug was causing the next line not to error
            f = lambda x: self.nll(x,h1.h_res.x[1],h2.h_res.x[1],data['Z1'],
                                   data['Z2'],data['score1'],data['score2'],
                                   data['scoreX'],data['N1'],data['N2'],
                                   self.M)
        except KeyError:
            f = lambda x: self.nll(x,h1.h_res.x[1],h2.h_res.x[1],data['Z1'],
                                   data['Z2'],data['score'],data['score'],
                                   data['score'],data['N1'],data['N2'],
                                   data.shape[0])
        try:
            data1['beta'] = data['beta1']
            data2['beta'] = data['beta2']
        except KeyError:
            pass
        pg = optimize.minimize_scalar(f,bounds=(-1.0,1.0),method='bounded',
                                     options={'disp':args.v})
        if not pg.success:
            sys.stderr.write(pg.message+'\n')
            raise ValueError
        if h1.h_res.x[1] == 0 or h2.h_res.x[1] == 0:
            pg.x = 0.0
        sy1 = self.estimate_sy(data1,h1.h_res.x[1],data1.shape[0])
        sy2 = self.estimate_sy(data2,h2.h_res.x[1],data2.shape[0])
        return h1, sy1, h2, sy2, pg, f

    def nll(self,pg,h1,h2,Z1,Z2,L1,L2,LX,N1,N2,M,W=None):
        Z = np.vstack((Z1,Z2)).T
        Z = Z.reshape((len(Z),2,1))
        S11 = 1 + (N1/M)*h1*L1
        S22 = 1 + (N2/M)*h2*L2
        sg = pg*np.sqrt(h1*h2)
        S12 = (np.sqrt(N1*N2)/M)*sg*LX
        S = np.dstack((np.vstack((S11,S12)).T,np.vstack((S12,S22)).T))
        logD = np.linalg.slogdet(S)[1]
        SIZ = np.linalg.solve(S,Z)
        ZTSIZ = np.sum(Z*SIZ,axis=(1,2))
        if W is None:
            l = np.sum(-0.5*logD - 0.5*ZTSIZ) - M*np.log(2*np.pi)
        else:
            l = np.sum(W*(-0.5*logD - 0.5*ZTSIZ)) - M*np.log(2*np.pi)
        return -1.0*l

    def convert_to_liability(self,K1,P1,K2,P2):
        self.res['Lia']=np.nan
        hl1 = self.h1_res.convert_to_liability(K1,P1)
        hl2 = self.h2_res.convert_to_liability(K2,P2)
        hgo=self.res['Val']['hg']
        tau1=stats.norm.ppf(1-K1)
        tau2=stats.norm.ppf(1-K2)
        hgl=hgo*(K1*(1-K1)*K2*(1-K2))/\
            (stats.norm.pdf(tau1)*stats.norm.pdf(tau2)*np.sqrt(P1*(1-P1)*P2*(1-P2)))
        self.res['Lia']['h1']=hl1
        self.res['Lia']['h2']=hl2
        self.res['Lia']['hg']=hgl
        self.res['Lia']['pg']=hgl/np.sqrt(hl1*hl2)
        return hgl

class fit_pg_pe(fit_pg):
    def __init__(self,data,args,M=None):
        if M is None: self.M = data.shape[0]
        else: self.M = M
        self.h1_res, self.sy1, self.h2_res, self.sy2, self.pg_res, self.ll =\
            self.__call__(data,args)
        res = np.array([self.h1_res.h_res.x[1], self.sy1, self.h2_res.h_res.x[1],
                        self.sy2, self.pg_res.x[0], self.pg_res.x[1]])
        self.null_ll = np.array([self.h1_res.null_ll[0], 0.0, self.h2_res.null_ll[0],
                                 0.0, self.ll((0.0,res[5])), self.ll((res[4],0.0))])
        self.alt_ll = np.array([self.h1_res.alt_ll[0], 0.0, self.h2_res.alt_ll[0], 0.0,
                                self.pg_res.fun, self.pg_res.fun])
        self.res = pd.DataFrame(np.vstack((res,np.tile(np.nan,len(res)))).T,
                                index=['h1','sy1','h2','sy2','pg','pe'],
                                columns=['Val','SE'])
        def close_call(x): # GET IT!!??!!
            res = self.__call__(x,args)
            return np.arr(res[0].h_res.x[1], res[1], res[2].h_res.x[1], res[3],
                    res[4].x[0], res[4].x[1])
        if not args.no_jackknife:
            self.jackknife = jackknife.jackknife(close_call, data, res, args)
            self.res['SE'] = self.jackknife.SE
        else:
            self.jackknife = None
        self.res['LR'] = 2*(self.null_ll-self.alt_ll)
        self.res['Z'] = self.res['Val']/self.res['SE']
        self.res['P (LRT)'] = 1-stats.chi2.cdf(self.res['LR'],1)
        self.res['P (Z)'] = 1-stats.chi2.cdf(self.res['Z']**2,1)

    def __call__(self,data,args):
        data1 = data.copy()
        data2 = data.copy()
        data1[['N','Z']] = data[['N1','Z1']]
        data2[['N','Z']] = data[['N2','Z2']]
        try:
            data1['beta'] = data['beta1']
            data2['beta'] = data['beta2']
        except KeyError:
            pass
        h1 = fit_h1(data1,args,jk=False,M=self.M)
        h2 = fit_h1(data2,args,jk=False,M=self.M)
        f = lambda x: self.nll(x,h1.h_res.x[1],h2.h_res.x[1],data['Z1'],data['Z2'],
                               data['score'],data['N1'],data['N2'],data['Ns'],
                               self.M)
        x0 = [0.01, 0.01] #self.initial_guess(f)
        pg = optimize.minimize(f,x0,bounds=((-1.0,1.0),(-1.0,1.0)),
                                     options={'disp':args.v})
        if not pg.success:
            sys.stderr.write(pg.message+'\n')
            raise ValueError
        sy1 = self.estimate_sy(data1,h1.h_res.x[1],data1.shape[0])
        sy2 = self.estimate_sy(data2,h2.h_res.x[1],data2.shape[0])
        return h1, sy1, h2, sy2, pg, f

    def initial_guess(self,f,start,finish,step,jac=True):
        # if jac: f0 = lambda x: f(x)[0]
        # else: f0 = f
        # p_ar = np.arange(start,finish+step,step)
        # f_ar = np.array(map(f0,p_ar))
        # index = np.argmin(f_ar)
        # return p_ar[index], f_ar[index]
        return 0.01,0.01

    def nll(self,x,h1,h2,Z1,Z2,L,N1,N2,Ns,M):
        pg,pe = x
        hg = pg*np.sqrt(h1*h2)
        he = pe*np.sqrt((1-h1)*(1-h2))
        Z = np.vstack((Z1,Z2)).T
        Z = Z.reshape((len(Z),2,1))
        S11 = 1 + (N1/M)*h1*L
        S22 = 1 + (N2/M)*h2*L
        S12 = (Ns/np.sqrt(N1*N2))*(he+hg) + (np.sqrt(N1*N2)/M)*hg*L
        S = np.dstack((np.vstack((S11,S12)).T,np.vstack((S12,S22)).T))
        logD = np.linalg.slogdet(S)[1]
        SIZ = np.linalg.solve(S,Z)
        ZTSIZ = np.sum(Z*SIZ,axis=(1,2))
        l = np.sum(-0.5*logD - 0.5*ZTSIZ) - M*np.log(2*np.pi)
        return -1.0*l
