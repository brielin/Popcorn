Popcorn
======

Popcorn is a program for estimaing the correlation of causal variant effect
sizes across populations in GWAS. Estimation involves two steps:

1. Compute mode: computing cross-population scores from a reference panel,
these roughly correspond to the similarity of the LD at each SNP across the populations
2. Fit mode: fitting the heritability and transethnic genetic correlation of
a pair of summary statistics files to the scores

There is no need to install Popcorn. Simply clone the git repository to a folder
called popcorn and run with python -m popcorn {compute,fit} <arguments> outfile. Please see the examples
below and the argument documentations for details.

Example usage:  
`python -m popcorn compute -v 1 --bfile1 /path/to/EUR_refpanel --bfile2 /path/to/EAS_refpanel scores.txt`  
`python -m popcorn fit -v 1 --cfile scores.txt --sfile1 /path/to/EUR_sumstats.txt --sfile2 /path/to/EAS_sumstats.txt EUR_EAS_corr.txt`  

For a full list of arguments and documentations type:  
`python -m popcorn compute -h`  
`python -m popcorn fit -h`  

Output:  
Popcorn reports the common-SNP observed scale (Val (Obs)) heritability (h^2) for both populations (h1^2 and h2^2),
the genetic effect or genetic impact correlation (pge or pgi), and the standard error of
these estimates (SE). If a transormation to the underlying liability scale is requested,
it also outputs the value on the liability scale (Val (Lia)). Popcorn also computes Z-scores
and p-values corresponding to these Z-scores. In the case of heritability, the p-value
reported is for a test that the heritability is *greater than 0*: P(h^2 > 0.0). In the
case of genetic correlation the p-value reported is for a test that the genetic correlation
is *less than 1.0*: P(pg < 1.0).

Dependences:  
Popcorn was developed using the following external python libraries.
If you experience problems running Popcorn, try updating your libraries,
epecially less stable ones like Pandas, Pysnptools and Statsmodels,
to the versions listed here.  
numpy 1.9.2  
scipy 0.16.1  
statsmodels 0.6.1  
pandas 0.17.1  
pysnptools 0.3.9  
bottleneck 1.0.0  

Test files:   
Popcorn comes with several files in the test directory to verify output. If you want
to make sure you have the software working correctly, run   
`python -m Popcorn compute -v 2 --bfile1 Popcorn/test/EUR_ref --bfile2 Popcorn/test/EAS_ref --gen_effect Popcorn/test/EUR_EAS_test_ge.cscore`   
The screen should show   
~~~
Popcorn version 0.9.6
(C) 2015-2016 Brielin C Brown
University of California, Berkeley
GNU General Public License v3

Invoking command: python /home/brielin/CoHeritability/Popcorn/__main__.py compute -v 2 --bfile1 Popcorn/test/EUR_ref --bfile2 Popcorn/test/EAS_ref Popcorn/test/EUR_EAS_test.cscore
Beginning analysis at DATE
50000 SNPs in file 1
50000 SNPs in file 2
49855 SNPs in file 1 after MAF filter
44021 SNPs in file 2 after MAF filter
43898 SNPs common in both populations
T1
T2
TC
Analysis finished at DATE
~~~

Then run   
`python -m Popcorn fit -v 1 --cfile Popcorn/test/EUR_EAS_test_ge.cscore --gen_effect --sfile1 Popcorn/test/EUR_test.txt --sfile2 Popcorn/test/EAS_test.txt Popcorn/test/EAS_EUR_test_ge_result.txt`   
The screen should show
~~~   
Popcorn version 0.9.6
(C) 2015-2016 Brielin C Brown
University of California, Berkeley
GNU General Public License v3

Invoking command: python /home/brielin/CoHeritability/Popcorn/__main__.py fit -v 1 --cfile Popcorn/test/EUR_EAS_test_ge.cscore --gen_effect --sfile1 Popcorn/test/EUR_test.txt --sfile2 Popcorn/test/EAS_test.txt Popcorn/test/EAS_EUR_test_ge_result.txt
Beginning analysis at DATE
Analyzing a pair of traits in two populations
49999 49999 SNPs detected in input.
43897 SNPs remaining after merging with score file.
43897 SNPs remaining after filtering on AF and self-compliment SNPs.
Initial estimate:
       Val (obs)  SE
h1^2   0.682685 NaN
h2^2   0.340257 NaN
pg     0.053720 NaN
Performing jackknife estimation of the standard error using 200 blocks of size 219 .
Jackknife iter: 200
      Val (obs)        SE          Z     P (Z)
h1^2   0.682685  0.050053  13.639243  0.000000
h2^2   0.340257  0.034388   9.894545  0.000000
pge    0.053720  3.836284   0.246666  0.805167
Analysis finished at DATE
Total time elapsed: T
~~~

Do not hesitate to contact me with questions as this is early stage software!

For more details of the method please see: http://biorxiv.org/content/early/2016/02/23/036657 or http://www.cell.com/ajhg/abstract/S0002-9297(16)30135-5