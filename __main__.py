#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from popcorn import fit
from popcorn import compute
from popcorn import sumstats
import numpy as np
import pandas as pd
import time
import sys
import argparse
import logging
import traceback
from IPython import embed


class Logger(object):
    def __init__(self,fname):
        self.terminal = sys.stdout
        self.log = open(fname+".log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

__version__='0.9.4'
header='Popcorn version '+__version__+'\n'\
'(C) 2015-2016 Brielin C Brown\n'\
'University of California, Berkeley\n'\
'GNU General Public License v3\n'

def main(args):
    # Arguments common to compute and fit modes
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('out',help='Specify output file.',
                               default='popcorn.out')
    parent_parser.add_argument('-v',help='Specify output verbosity.',type=int,
                               choices=[0,1,2,3],default=0)
    parent_parser.add_argument('--maf',help='Specify MAF cutoff. We suggest'
                               ' at least 0.01.',type=float, default=0.05)
    parent_parser.add_argument('--gen_effect',help='Compute genetic effect'
                               ' correlation instead of genetic impact'
                               ' correlation. This must be specified both'
                               ' when computing scores and fitting the model.',
                               default=False,action='store_true')
    parent_parser.add_argument('--from_bp',default=None,type=int,help='Specify'
                               ' base position to start analysis.')
    parent_parser.add_argument('--to_bp',default=None,type=int,help='Specify'
                               ' base position to finish analysus.')
    parent_parser.add_argument('--no_align',default=False,action='store_true',
                               help=argparse.SUPPRESS)

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='Program mode: "compute" for'
                                       ' computing covariance scores, "fit" for'
                                       ' fitting summary statistics with pre-'
                                       'computed covariance scores.',
                                       dest='mode')
    # Arguments exclusive to compute mode
    parser_compute = subparsers.add_parser('compute',parents=[parent_parser])
    parser_compute.add_argument('--bfile',help='Specify binary plink filename,'
                                ' without extension',
                                default=None)
    parser_compute.add_argument('--bfile1',help='Specify binary plink filename '
                                'for population 1, without extension',default=None)
    parser_compute.add_argument('--bfile2',help='Specify binary plink filename '
                                'for population 2, without extension',default=None)
    parser_compute.add_argument('--window_size',help='Specify window size,'
                                'in units of WINDOW_TYPE',
                                default=1000,type=int)
    parser_compute.add_argument('--window_type',help='Specify window type'
                                ' (SNP,KBP)',default='KBP')
    parser_compute.add_argument('--SNPs_to_read',help='Specify number of SNPs'
                                ' to read at a time. Mainly for debugging.',
                                default=1000,type=int)
    parser_compute.add_argument('--SNPs_to_store',help='Specify size of'
                                ' in memory SNP array. May need to increase'
                                ' for very dense panels or wide window sizes.',
                                type=int,default=10000)
    parser_compute.add_argument('--extract',default=None,help='Specify a text'
                                ' file containing a whitespace separated list'
                                ' of SNP IDs for use in the computation.')
    # parser_compute.add_argument('--afile',help='Specify ancestry file for'
    #                             'conditioning in admixed populations.',
    #                             default=None)
    # parser_compute.add_argument('--afile1',help='Specify ancestry file for'
    #                             'conditioning in admixed populations.',
    #                             default=None)
    # parser_compute.add_argument('--afile2',help='Specify ancestry file for'
    #                             'conditioning in admixed populations.',
    #                             default=None)
    # Arguments exclusive to fit mode
    parser_fit = subparsers.add_parser('fit',parents=[parent_parser])
    parser_fit.add_argument('--sfile',help='Specify tab delimited'
                            ' summary statistics file')
    parser_fit.add_argument('--sfile1',help='Specify tab delimited'
                            ' summary statistics file for population 1')
    parser_fit.add_argument('--sfile2',help='Specify tab delimited'
                            ' summary statistics file for population 2')
    parser_fit.add_argument('--cfile',help='Specify covariance scores file')
    parser_fit.add_argument('--Ns',help='Either 1) the number'
                            ' of overlapping samples or 2) proportion of'
                            ' overlapping samples or a 3) text file containing'
                            ' the number of overlapping samples per SNP ID.')
    parser_fit.add_argument('--regions',help=argparse.SUPPRESS)
    parser_fit.add_argument('--tol',help='Specify the convergence tolerance'
                            ' for MLE. Unless you have convergence issues use'
                            ' default',default=.00001)
    parser_fit.add_argument('--no_jackknife',default=False,action='store_true',
                            help='Just compute the point estimate without computing'
                            ' the jackknife standard error.')
    parser_fit.add_argument('--K1',default=None,type=float,help='Specify'
                            ' population prevelance of binary phenotype 1.')
    parser_fit.add_argument('--P1',default=None,type=float,help='Specify'
                            ' study prevelance of binary phenotype 1.')
    parser_fit.add_argument('--K2',default=None,type=float,help='Specify'
                            ' population prevelance of binary phenotype 2.')
    parser_fit.add_argument('--P2',default=None,type=float,help='Specify'
                            ' study prevelance of binary phenotype 2.')
    parser_fit.add_argument('--M',default=None,type=int,help=argparse.SUPPRESS)
    parser_fit.add_argument('--no_intercept',default=False,action='store_true',
                            help='Fixes intercept of heritability to be 1.0.')
    args = parser.parse_args()

    # Set up logger to print to log file and std out
    # logging.basicConfig(level=logging.INFO, format='%(message)s')
    # logger = logging.getLogger()
    # logger.addHandler(logging.FileHandler(args.out+'.log'))
    # print = logger.info

    sys.stdout = Logger(args.out+'.o')
    sys.stderr = Logger(args.out+'.e')
    print(header)
    print('Invoking command: python '+' '.join(sys.argv))
    print('Beginning analysis at {T}'.format(T=time.ctime()))
    start_time = time.time()
    try:
        if args.mode == 'fit':
            if args.cfile is None:
                raise ValueError('Must provide summary covariance scores.')
            if (args.K1 is None) ^ (args.P1 is None):
                raise ValueError('Must provide K and P to convert'
                                 'observed scale to liability scale.')
            if (args.sfile is not None) \
                    ^ (args.sfile1 is not None) \
                    and (args.sfile2 is None):
                print("Analyzing a single trait in one population")
                if args.sfile is None: args.sfile = args.sfile1
                # with open(args.cfile,'r') as f:
                #     M = int(f.readline().strip().split()[-1])
                scores = pd.read_table(args.cfile,header=None,comment='#',
                                       names=['chr','pos','id','a1','a2','af',
                                               'score'])
                M = scores.shape[0]
                scores.index = scores['id']
                data = sumstats.sumstats_1_trait(scores,args)
                if (args.M is not None) and (args.M<M):
                    data.data = data.data.iloc[range(args.M)]
                if args.regions:
                    res = fit.fit_by_region(data.data,args,t='h1')
                else:
                    res = fit.fit_h(data.data,args,M=M)
            elif (args.sfile is not None) \
                    ^ (args.sfile1 is not None) \
                    and (args.sfile2 is not None):
                if args.sfile1 is None: args.sfile1 = args.sfile
                with open(args.cfile,'r') as f:
                    for line in f:
                        if line[0]!='#':
                            ncols = len(line.strip().split())
                            break
                if ncols == 7:
                    print("Analyzing a pair of traits in one population")
                    scores = pd.read_table(args.cfile,header=None,comment='#',
                                           names=['chr','pos','id','a1','a2',
                                                  'af','score'],sep='\t')
                elif ncols == 10:
                    print("Analyzing a pair of traits in two populations")
                    scores = pd.read_table(args.cfile,header=None,comment='#',
                                           names=['chr','pos','id','a1','a2',
                                                  'af1','af2','score1','score2',
                                                  'scoreX'])
                else:
                    raise ValueError('Cscores file not understood. Did you'
                                     ' compute them with this program?')
                M = scores.shape[0]
                scores.index = scores['id']
                data = sumstats.sumstats_2_trait(scores,args)
                if (args.M is not None) and (args.M<M):
                    data.data = data.data.iloc[range(args.M)]
                if data.overlap:
                    if args.regions:
                        res = fit.fit_by_region(data.data,args,t='pg_pe',M=M)
                    else:
                        res = fit.fit_pg_pe(data.data,args,M=M)
                else:
                    if args.regions:
                        res = fit.fit_by_region(data.data,args,t='pg',M=M)
                    else:
                        res = fit.fit_pg(data.data,args,M=M)
            else:
                raise ValueError('Must provide sfile or sfile1, or sfile1 and 2')
            res.write(args.out)
        elif args.mode == 'compute':
            if (args.bfile is not None) \
                    ^ (args.bfile1 is not None) \
                    and (args.bfile2 is None):
                if args.bfile is None: args.bfile = args.bfile1
                scores = compute.covariance_scores_1_pop(args)
            elif (args.bfile is not None) \
                    ^ (args.bfile1 is not None) \
                    and (args.bfile2 is not None):
                if args.bfile1 is None: args.bfile1 = args.bfile
                scores = compute.covariance_scores_2_pop(args)
            else:
                raise ValueError('Must provide bfile or bfile1, or bfile1 and 2')
            scores.write(args)
    except Exception:
        # ex_type, ex, tb = sys.exc_info()
        # traceback.print_last()
        # print(traceback.format_exc(ex))
        raise
    print('Analysis finished at {T}'.format(T=time.ctime()) )
    time_elapsed = round(time.time()-start_time,2)
    print('Total time elapsed: {T}'.format(T=time_elapsed))

if __name__ == '__main__':
    main(sys.argv)
