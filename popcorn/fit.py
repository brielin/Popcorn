from __future__ import division
from __future__ import print_function
import numpy as np
import pandas as pd
import sys
import argparse
import logging
import traceback
import statsmodels.api as sm
from scipy import optimize, stats
# from IPython import embed
from time import time
from collections import namedtuple
from popcorn import jackknife

class fit_by_region(object):
    def __init__(self,data,args,t,M):
        types = {'h1':fit_h,'pg':fit_pg,'pg_pe':fit_pg_pe}
        regions = pd.read_table(args.regions,sep='\s*')
        if regions['chr'][0][0]=='c':
            regions['chr'] = map(lambda x: int(x[3:]), regions['chr'])
        regions['M'] = self.get_M_by_region(regions,data)
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
                #self.__dict__[k].to_csv(outfile+'.'+k,sep='\t',na_rep='NaN')
                f=open(outfile+'.'+k,'w')
                print(self.__dict__[k].to_string())
                print(self.__dict__[k].to_string(),file=f)
            except AttributeError:
                pass


class fit_h(object):
    def __init__(self,data,args,K=None,P=None,jk=True,M=None):
        if M is None: self.M = data.shape[0]
        else: self.M = M
        self.h_res, self.sy, self.ll = self.__call__(data,args)
        #print(self.h_res.x)
        # if args.no_intercept:
        #     self.null_ll = np.array([self.ll(0.0), 0.0])
        # else:
        #     # This is still wrong though
        #     self.null_ll = np.array([self.ll([1.0,0.0]), 0.0])
        # self.alt_ll = np.array([self.h_res.fun, 0.0])
        res = np.array([self.h_res.x[1], self.sy])
        #res = np.array([self.h_res.x, self.sy])
        self.res = pd.DataFrame(np.vstack((res,np.tile(np.nan,len(res)))).T,
                                index=['h^2','sy'],columns=['Val (obs)','SE'])
        def close_call(x):
            res = self.__call__(x,args)
            return res[0].x[1], res[1]
            #return res[0].x, res[1]
        if jk and not args.no_jackknife:
            self.jackknife = jackknife.jackknife(close_call, data, res, args)
            self.res['SE'] = self.jackknife.SE
        else:
            self.jackknife = None
        if (K is not None)&(P is not None):
            self.convert_to_liability(K,P)
        self.res['Z'] = self.res['Val (obs)']/self.res['SE']
        self.res['P(h^2 > 0)'] = 1-stats.chi2.cdf(self.res['Z']**2,1)
        self.res = pd.DataFrame(self.res.loc['h^2'])
        # if args.plot_likelihood:
        #     if args.no_intercept:
        #         f=lambda x: self.nll_no_intercept(x,data['Z'],data['N'],self.M,
        #                                           data['score'],W=W)
        #     else:
        #         f=lambda x: self.nll(x,data['Z'],data['N'],self.M,data['score'],W=W)
        #     x_vals = np.linspace(0, 1, num=100)
        #     l_vals = np.array([f(x) for x in x_vals])
        #     plt.plot(x_vals, l_vals, 'bs--')
        # self.res['LR'] = 2*(self.null_ll-self.alt_ll)
        # self.res['P (LRT)'] = 1-stats.chi2.cdf(self.res['LR'],1)

    def __call__(self,data,args):
        W = 1.0/np.maximum(data['score'],np.ones(data['score'].shape))
        if args.no_intercept:
            if not args.use_mle:
                f = None
                Y = (data['Z']**2)
                X = data['score']*(data['N']/self.M)
                #beta, _, _, _ = np.linalg.lstsq(X.reshape((X.shape[0], 1)), Y - 1.)
                wls = sm.WLS(Y, X, W)
                res_wls = wls.fit()
                H_Res = namedtuple('H_Res', ['x'])
                h = H_Res(x=[1.0, res_wls.params.values[0]])
                sy = self.estimate_sy(data,h.x[1],data.shape[0])
            else:
                f=lambda x: self.nll_no_intercept(x,data['Z'],data['N'],self.M,
                                                  data['score'],W=W)
                h = optimize.minimize_scalar(f,bounds=(0.0,1.0),method='bounded',
                                             options={'disp':args.v-1,'xatol':args.tol})
                sy = self.estimate_sy(data,h.x,data.shape[0])
                h.x = [1.0,h.x]
        else:
            if not args.use_mle:
                f = None
                Y = (data['Z']**2)
                X = data['score']*(data['N']/self.M)
                X = np.stack([np.ones(X.shape), X.values], 1)
                #beta, _, _, _ = np.linalg.lstsq(X, Y)
                wls = sm.WLS(Y, X, W)
                res_wls = wls.fit()
                H_Res = namedtuple('H_Res', ['x'])
                h = H_Res(x=list(res_wls.params))
                sy = self.estimate_sy(data,h.x[1],data.shape[0])
            else:
                f=lambda x: self.nll(x,data['Z'],data['N'],self.M,data['score'],W=W)
                for x1 in np.arange(0.0,1,0.05):
                    x0=[1.0,x1]
                    h = optimize.minimize(f,x0,bounds=((None,None),(0.0,1.0)),
                                          options={'disp':args.v-1})
                    if not h.success:
                        continue
                    else:
                        break
                if not h.success:
                    sys.stderr.write(h.message+'\n')
                    h.x = [np.nan,np.nan]
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

    def nll_no_intercept(self,x,Z,N,M,score,W=None):
        S = (1 + (N/M)*x*score)
        if W is None:
            l = -0.5*M*np.log(2*np.pi) - 0.5*np.sum(np.log(S))\
                - 0.5*np.sum(Z**2/S)
        else:
            l = -0.5*M*np.log(2*np.pi) - 0.5*np.sum(W*np.log(S))\
                - 0.5*np.sum(W*(Z**2/S))
        return -1.0*l

    def write(self,outfile):
        print(self.res.to_string())
        self.res.to_csv(outfile,sep='\t',na_rep='NaN')

    def convert_to_liability(self,K,P):
        self.res['Val (Lia)']=np.nan
        ho=self.res['Val (obs)']['h^2']
        tau=stats.norm.ppf(1-K)
        hl=ho*((K*(1-K))**2)/(stats.norm.pdf(tau)**2*P*(1-P))
        self.res['Val (Lia)']['h^2']=hl
        return hl

class fit_pg(fit_h):
    def __init__(self,data,args,M=None):
        if M is None: self.M = data.shape[0]
        else: self.M = M
        self.h1_res, self.sy1, self.h2_res, self.sy2, self.pg_res, self.ll =\
            self.__call__(data,args)
        if args.plot_likelihood:
            # This is bad, but it will avoid the dependency for people that dont
            #  use this flag, most people shouldn't.
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            x_vals = np.linspace(-1.5, 1.5, num=200)
            l_vals = np.array([self.ll(x) for x in x_vals])
            plt.plot(x_vals, l_vals, 'bs--')
            plt.savefig(args.out + '.png')

        h1r = self.h1_res.h_res.x[1]
        h2r = self.h2_res.h_res.x[1]
        pgr = self.pg_res.x
        hgr = pgr*np.sqrt(h1r*h2r)
        res = np.array([h1r, self.sy1, h2r, self.sy2, hgr, pgr])
        # self.null_ll = np.array([self.h1_res.null_ll[0], 0.0,
        #                          self.h2_res.null_ll[0], 0.0,
        #                          self.ll(1.0),self.ll(0.0)])
        # self.alt_ll = np.array([self.h1_res.alt_ll[0], 0.0,
        #                         self.h2_res.alt_ll[0], 0.0,
        #                         self.pg_res.fun,self.pg_res.fun])
        self.res = pd.DataFrame(np.vstack((res,np.tile(np.nan,len(res)))).T,
                                index=['h1^2','sy1','h2^2','sy2','hg','pg'],
                                columns=['Val (obs)','SE'])
        if args.v > 0: print("Initial estimate:\n",self.res.loc[['h1^2','h2^2','pg']])
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
        self.res['Z'] = (self.res['Val (obs)'])/self.res['SE']
        self.res['Z']['pg'] = (1.0-self.res['Val (obs)']['pg'])/self.res['SE']['pg']
        self.res['P (Z)'] = 1-stats.chi2.cdf(self.res['Z']**2,1)
        if args.gen_effect:
            self.res.index=['h1^2','sy1','h2^2','sy2','hge','pge']
            self.res=self.res.loc[['h1^2','h2^2','pge']]
        else:
            self.res.index=['h1^2','sy1','h2^2','sy2','hgi','pgi']
            self.res=self.res.loc[['h1^2','h2^2','pgi']]

    def __call__(self,data,args):
        t=time()
        data1 = data.copy()
        data2 = data.copy()
        try: #two populations
            data1[['N','af','Z','score']] = data[['N1','af1','Z1','score1']]
            data2[['N','af','Z','score']] = data[['N2','af2','Z2','score2']]
            W = 1.0/np.maximum(data['scoreX'],np.ones(data['scoreX'].shape))
        except KeyError: #one population
            data1[['N','Z']] = data[['N1','Z1']]
            data2[['N','Z']] = data[['N2','Z2']]
            W = 1.0/np.maximum(data['score'],np.ones(data['score'].shape))
        h1 = fit_h(data1,args,args.K1,args.P1,jk=False,M=self.M)
        h2 = fit_h(data2,args,args.K2,args.P2,jk=False,M=self.M)
        if 'score1' in data.columns:
            L1=data['score1']
            L2=data['score2']
            LX=data['scoreX']
        else:
            L1=data['score']
            L2=data['score']
            LX=data['score']
        if not args.use_mle:
            Y = data['Z1']*data['Z2']
            X = LX*np.sqrt(data['N1']*data['N2'])/(self.M)
            if not args.no_intercept:
                X = sm.add_constant(X)
            wls = sm.WLS(Y, X, W)
            res_wls = wls.fit()
            P_Res = namedtuple('P_Res', ['x'])
            pg = P_Res(x=res_wls.params.values[1]/np.sqrt(h1.h_res.x[1]*h2.h_res.x[1]))
            f=None
        else:
            f = lambda x: self.nll(x,h1.h_res.x[1],h2.h_res.x[1],data['Z1'],
                                   data['Z2'],L1,L2,LX,data['N1'], data['N2'],
                                   self.M,W=W)
            try:
                data1['beta'] = data['beta1']
                data2['beta'] = data['beta2']
            except KeyError:
                pass
            pg = optimize.minimize_scalar(f,bounds=(-1.0,1.0),method='bounded',
                                        options={'disp':args.v-1,'xatol':args.tol})
            if not pg.success:
                sys.stderr.write(pg.message+'\n')
                pg.x=np.nan

        if h1.h_res.x[1] == 0 or h2.h_res.x[1] == 0:
            pg.x = np.nan
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
        self.res['Val (Lia)']=np.nan
        hl1 = self.h1_res.res['h^2']['Val (Lia)']
        hl2 = self.h2_res.res['h^2']['Val (Lia)']
        hgo=self.res['Val (obs)']['hg']
        tau1=stats.norm.ppf(1-K1)
        tau2=stats.norm.ppf(1-K2)
        hgl=hgo*(K1*(1-K1)*K2*(1-K2))/\
            (stats.norm.pdf(tau1)*stats.norm.pdf(tau2)*np.sqrt(P1*(1-P1)*P2*(1-P2)))
        self.res['Val (Lia)']['h1^2']=hl1
        self.res['Val (Lia)']['h2^2']=hl2
        self.res['Val (Lia)']['hg']=hgl
        self.res['Val (Lia)']['pg']=hgl/np.sqrt(hl1*hl2)
        return hgl

class fit_pg_pe(fit_pg):
    def __init__(self,data,args,M=None):
        if M is None: self.M = data.shape[0]
        else: self.M = M
        self.h1_res, self.sy1, self.h2_res, self.sy2, self.pg_res, self.ll =\
            self.__call__(data,args)
        res = np.array([self.h1_res.h_res.x[1], self.sy1, self.h2_res.h_res.x[1],
                        self.sy2, self.pg_res.x[0], self.pg_res.x[1]])
        # self.null_ll = np.array([self.h1_res.null_ll[0], 0.0, self.h2_res.null_ll[0],
        #                          0.0, self.ll((0.0,res[5])), self.ll((res[4],0.0))])
        # self.alt_ll = np.array([self.h1_res.alt_ll[0], 0.0, self.h2_res.alt_ll[0], 0.0,
        #                         self.pg_res.fun, self.pg_res.fun])
        self.res = pd.DataFrame(np.vstack((res,np.tile(np.nan,len(res)))).T,
                                index=['h1','sy1','h2','sy2','pg','pe'],
                                columns=['Val (obs)','SE'])
        def close_call(x): # GET IT!!??!!
            res = self.__call__(x,args)
            return np.arr(res[0].h_res.x[1], res[1], res[2].h_res.x[1], res[3],
                    res[4].x[0], res[4].x[1])
        if not args.no_jackknife:
            self.jackknife = jackknife.jackknife(close_call, data, res, args)
            self.res['SE'] = self.jackknife.SE
        else:
            self.jackknife = None
        self.res['Z'] = self.res['Val (obs)']/self.res['SE']
        self.res['P (Z)'] = 1-stats.chi2.cdf(self.res['Z']**2,1)
        if args.gen_effect:
            self.res.index=['h1^2','sy1','h2^2','sy2','hge','pge']
            self.res=self.res.loc[['h1^2','h2^2','pge']]
        else:
            self.res.index=['h1^2','sy1','h2^2','sy2','hge','pgi']
            self.res=self.res.loc[['h1^2','h2^2','pgi']]
        # self.res['LR'] = 2*(self.null_ll-self.alt_ll)
        # self.res['P (LRT)'] = 1-stats.chi2.cdf(self.res['LR'],1)

    def __call__(self,data,args):
        data1 = data.copy()
        data2 = data.copy()
        data1[['N','Z']] = data[['N1','Z1']]
        data2[['N','Z']] = data[['N2','Z2']]
        W = 1.0/np.maximum(data['score'],np.ones(data['score'].shape))
        try:
            data1['beta'] = data['beta1']
            data2['beta'] = data['beta2']
        except KeyError:
            pass
        h1 = fit_h(data1,args,jk=False,M=self.M)
        h2 = fit_h(data2,args,jk=False,M=self.M)
        f = lambda x: self.nll(x,h1.h_res.x[1],h2.h_res.x[1],data['Z1'],data['Z2'],
                               data['score'],data['N1'],data['N2'],data['Ns'],
                               self.M,W=W)
        x0 = [0.01, 0.01] #self.initial_guess(f)
        pg = optimize.minimize(f,x0,bounds=((-1.0,1.0),(-1.0,1.0)),
                                     options={'disp':args.v-1})
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
