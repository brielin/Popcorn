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

class fit_h1(object):
    def __init__(self,data,args):
        #self.h1 = self.newtons_method(data,args,data.shape[0],0.0,1.0,0.05)
        f=lambda x: self.nll(x,data['Z'],data['N'],data.shape[0],data['score'])
        # h0,l0 = self.initial_guess(f,0.0,1.0,0.05,jac=True)
        # h = optimize.minimize(f,h0,jac=True,bounds=(0.0,1.0),
        #                       options={'disp':True})
        h = optimize.minimize_scalar(f,bounds=(0.0,1.0),method='bounded',
                              options={'disp':True})
        self.h_res = h
        self.sy = self.estimate_sy(data,h.x,data.shape[0])
        # Jackknife

    def estimate_sy(self,data,h1,M):
        try:
            B = data['beta']
        except KeyError:
            return np.nan
        var = 2*data['af']*(1-data['af'])
        A = (1/(data['N']*var))*(1 + (data['N']/M)*h1*data['score'])
        return np.sum(B**2/A)/M

    def nll(self,x,Z,N,M,score):
        S = (1 + (N/M)*x*score)
        #dSdh1 = (N/M)*score
        #fp = np.sum(-dSdh1/S + (dSdh1/S**2)*Z**2)
        l = -0.5*M*np.log(2*np.pi) - 0.5*np.sum(np.log(S))\
            - 0.5*np.sum(Z**2/S)
        return -1.0*l
        #return (-1.0*l, -1.0*fp)

class fit_pg(fit_h1):
    def __init__(self,data,args):
        data1 = data.copy()
        data2 = data.copy()
        data1['N'] = data['N1']
        data2['N'] = data['N2']
        data1['af'] = data['af1']
        data2['af'] = data['af2']
        data1['Z'] = data['Z1']
        data2['Z'] = data['Z2']
        data1['score'] = data['score1']
        data2['score'] = data['score2']
        try:
            data1['beta'] = data['beta1']
            data2['beta'] = data['beta2']
        except KeyError:
            pass
        h1 = fit_h1(data1,args)
        h2 = fit_h1(data2,args)
        f = lambda x: self.nll(x,data,h1.h_res.x,h2.h_res.x,data.shape[0])
        pg = optimize.minimize_scalar(f,bounds=(0.0,1.0),method='bounded',
                                     options={'disp':True})
        self.pg_res = pg
        self.h1_res = h1
        self.h2_res = h2
        self.sy1 = self.estimate_sy(data,h1.h_res.x,data1.shape[0])
        self.sy2 = self.estimate_sy(data,h2.h_res.x,data2.shape[0])

    def nll(self,pg,data,h1,h2,M):
        Z = np.vstack((data['Z1'],data['Z2'])).T.reshape((M,2,1))
        S11 = 1 + (data['N1']/M)*h1*data['score1']
        S22 = 1 + (data['N2']/M)*h2*data['score2']
        sg = pg*np.sqrt(h1*h2)
        S12 = (np.sqrt(data['N1']*data['N2'])/M)*sg*data['scoreX']
        S = np.dstack((np.vstack((S11,S12)).T,np.vstack((S12,S22)).T))
        logD = np.linalg.slogdet(S)[1]
        SIZ = np.linalg.solve(S,Z)
        ZTSIZ = np.sum(Z*SIZ,axis=(1,2))
        l = np.sum(-0.5*logD - 0.5*ZTSIZ) - M*np.log(2*np.pi)
        return -1.0*l


    # def ll(pg,data,h1,h2,M):
    #     Z = np.vstack((data['Z1'],data['Z2'])).T.reshape((M,2,1))
    #     S11 = 1 + (data['N1']/M)*h1*data['score1']
    #     S22 = 1 + (data['N2']/M)*h2*data['score2']
    #     sg = pg*np.sqrt(h1*h2)
    #     S12 = (np.sqrt(score['N1']*score['N2'])/M)*sg*data['scoreX']
    #     S = np.dstack((np.vstack((S11,S12)).T,np.vstack((S12,S22)).T))
    #     dS12dpg = ((np.sqrt(score['N1']*score['N2'])/M)*np.sqrt(h1*h2)*\
    #                    score['scoreX']).reshape((M,1))
    #     dSdpg1 = np.hstack((np.zeros(dS12dpg.shape),dS12dpg)).reshape((M, 2,1))
    #     dSdpg = np.dstack((np.hstack((np.zeros(dS12dpg.shape),dS12dpg)),
    #                        np.hstack((dS12dpg,np.zeros(dS12dpg.shape)))))
    #     logD = np.linalg.slogdet(S)[1]
    #     SIZ = np.linalg.solve(S,Z)
    #     ZTSIZ = np.sum(Z*SIZ,axis=(1,2))
    #     SIdSdpg = np.linalg.solve(S,dSdpg1)
    #     TrSIdSdpg = 2*SIdSdpg[:,0,0]
    #     dSdpgSIZ = np.sum(dSdpg*SIZ,axis=1).reshape((M,2,1))
    #     ZTSIdSdpgSIZ = np.sum(SIZ*dSdpgSIZ,axis=(1,2))
    #     TrSIdSdpgdSdpg = 2*np.sum(SIdSdpg*SIdSdpg,axis=(1,2))
    #     l = np.sum(-0.5*logD - 0.5*ZTSIZ) - M*np.log(2*np.pi)
    #     g = np.sum(-0.5*TrSIdSdpg + 0.5*ZTSIdSdpgSIZ)
    #     I = np.sum(-0.5*TrSIdSdpgdSdpg)
    #     return g,I,l

    def initial_guess(self,f,start,finish,step,jac=True):
        if jac: f0 = lambda x: f(x)[0]
        else: f0 = f
        p_ar = np.arange(start,finish+step,step)
        f_ar = np.array(map(f0,p_ar))
        index = np.argmin(f_ar)
        return p_ar[index], f_ar[index]

class fit_pg_pe(fit_pg):
    def __init__(self,data,args):
        data1 = data.copy()
        data2 = data.copy()
        data1['N'] = data['N1']
        data2['N'] = data['N2']
        data1['af'] = data['af1']
        data2['af'] = data['af2']
        data1['Z'] = data['Z1']
        data2['Z'] = data['Z2']
        data1['score'] = data['score1']
        data2['score'] = data['score2']
        try:
            data1['beta'] = data['beta1']
            data2['beta'] = data['beta2']
        except KeyError:
            pass
        h1 = fit_h1(data1,args)
        h2 = fit_h1(data2,args)
        f = lambda x: self.nll(x,data,h1.h_res.x,h2.h_res.x,data.shape[0])
        pg = optimize.minimize_scalar(f,bounds=(0.0,1.0),method='bounded',
                                     options={'disp':True})
        self.pg_res = pg
        self.h1_res = h1
        self.h2_res = h2
        self.sy1 = self.estimate_sy(data,h1.h_res.x,data1.shape[0])
        self.sy2 = self.estimate_sy(data,h2.h_res.x,data2.shape[0])

    pass
