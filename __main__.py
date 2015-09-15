from __future__ import division
from __future__ import print_function
from my_tool import fit
from my_tool import compute
from my_tool import sumstats
import numpy as np
import pandas as pd
import time
import sys
import argparse
import logging
import traceback
from IPython import embed

def main(args):
    # Arguments common to compute and fit modes
    parser = argparse.ArgumentParser()
    parser.add_argument('-v',help='Specify output verbosity',type=int,
                        choices=[0,1,2,3])
    parser.add_argument('--out',help='Specify output file',default='my_tool.out')
    parser.add_argument('--maf',help='Specify MAF cutoff',type=float,
                        default=0.01)
    parser.add_argument('--per_allele',help='Use the per-allele effects model'
                        ' rather than the standardized effects model.',
                        default=False,action='store_true')
    subparsers = parser.add_subparsers(help='Program mode: "compute" for '
                                       'computing covariance scores, "fit" for '
                                       'fitting summary statistics with already'
                                       ' computed covariance scores.',
                                       dest='mode')
    # Arguments exclusive to compute mode
    parser_compute = subparsers.add_parser('compute')
    parser_compute.add_argument('--bfile',help='Speficy binary plink filename.',
                                default=None)
    parser_compute.add_argument('--bfile1',help='Specify binary plink filename'
                                'for population 1',default=None)
    parser_compute.add_argument('--bfile2',help='Specify binary plink filename'
                                'for population 2',default=None)
    parser_compute.add_argument('--window_size',help='Specify window size',
                                default=1000)
    parser_compute.add_argument('--window_type',help='Specify window type'
                                ' (SNP,BP)',default='SNP')
    parser_compute.add_argument('--SNPs_to_read',help='Specify number of SNPs'
                                'to read at a time. Do not set yourself',
                                default=1000)
    # Arguments exclusive to fit mode
    parser_fit = subparsers.add_parser('fit')
    parser_fit.add_argument('--sfile',help='Specify summary statistics file')
    parser_fit.add_argument('--sfile1',help='Specify summary statistics file')
    parser_fit.add_argument('--sfile2',help='Specify summary statistics file')
    parser_fit.add_argument('--cfile',help='Specify covariance scores file')
    parser_fit.add_argument('--Ns',help='Either 1) the number'
                            ' of overlapping samples or 2) proportion of'
                            ' overlapping samples or a 3) text file containing'
                            ' the number of overlapping samples per SNP ID.')
    parser_fit.add_argument('--regions',help='A text file with three columns'
                            ' specifying chromosome, starting base and ending'
                            ' base for a list of regions to analyse separately.')
    # parser_fit.add_argument('--niter',help='Specify the max iterations for MLE.'
    #                         ' Unless you have convergence issues use default',
    #                         default=1000,type=int)
    # parser_fit.add_argument('--ll_tol',help='Specify the convergence tolerance'
    #                         ' for MLE. Unless you have convergence issues use'
    #                         ' default',default=0.001)
    args = parser.parse_args()

    # Set up logger to print to log file and std out
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logger = logging.getLogger()
    logger.addHandler(logging.FileHandler(args.out+'.log'))
    print = logger.info

    print('Beginning analysis at {T}'.format(T=time.ctime()))
    start_time = time.time()
    try:
        if args.mode == 'fit':
            if args.sfile is None or args.cfile is None:
                raise ValueError('Must provide summary statistics and'
                                 ' covariance scores.')
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
                if args.regions:
                    res = fit.fit_by_region(data.data,scores,args,t='h1')
                else:
                    res = fit.fit_h1(data.data,args)
            elif (args.sfile is not None) \
                    ^ (args.sfile1 is not None) \
                    and (args.sfile2 is not None):
                if args.sfile1 is None: args.sfile1 = args.sfile
                with open(args.cfile,'r') as f:
                    ncols = len(f.readline().strip().split())
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
                raise ValueError('Must provide bfile or bfile1, or bfile1 and 2')
            #print(res.res.to_string())
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
        ex_type, ex, tb = sys.exc_info()
        print(traceback.format_exc(ex))
        raise
    print('Analysis finished at {T}'.format(T=time.ctime()) )
    time_elapsed = round(time.time()-start_time,2)
    print('Total time elapsed: {T}'.format(T=time_elapsed))

if __name__ == '__main__':
    main(sys.argv)
