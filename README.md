Popcorn
======

Popcorn is a program for estimaing the correlation of causal variant effect. This is the python3 version of Popcorn and still under development
sizes across populations in GWAS. Estimation involves two steps:

1. Compute mode: computing cross-population scores from a reference panel,
these roughly correspond to the similarity of the LD at each SNP across the populations
2. Fit mode: fitting the heritability and transethnic genetic correlation of
a pair of summary statistics files to the scores

# NOTE
Popcorn was written in Python 2.7 and Python 3+ is not supported. Please note this software was written in 2015, and I am not currently working on any projects using this tool for analysis. Thus, I do not anticipate additional development of this package in the near future. If you would like to or have ported this code to Python3, or you have generally found/fixed bugs, I am always happy to review a pull request. I do try to respond to emails or issues on here, but please note all maintence of this package happens in my spare time therefore I may be (sometimes very) slow to respond.



# Installation
Popcorn can be installed with setuptools or it can be run by executing the main. To install with setuptools

```
cd Popcorn
python setup.py install
```

and all of the requirements will be installed automatically. In this case you can then run popcorn
with 

`popcorn {compute, fit} <arguments> outfile`

Otherwise, install the requirements (using `pip install -r requirements.txt` or otherwise) and
execute the main directly with `python Popcorn/popcorn/__main__.py {compute, fit} <arguments> outfile`

Please see the examples below and the argument documentations for details.

# File formats and usage information
For computing scores, reference genotypes in plink binary (BED) format are required. See <http://www.cog-genomics.org/plink/1.9/input#bed> for more information on that format. PLINK2 can also be used for converting reference genomes in VCF format to BED format. The basic command to compute scores for two populations is

`popcorn compute -v 1 --bfile1 /path/to/EUR_refpanel --bfile2 /path/to/EAS_refpanel scores.txt`  

Which will create a file "scores.txt" in the local directory with 10 columns and no header. In order, the columns represent Chromosome, Base Position, SNP ID, Allele 1, Allele 2, Frequency of A1 in POP1, Frequency of A1 in POP2, LD score in POP1, LD score in POP2, and Cross-covariance score. If you compute scores for each chromosome separately, you can easily cat them together:

`cat scores_{1..22}.txt > scores_all.txt`

Pre-computed scores for EUR and EAS 1000 genomes populations are provided at <https://www.dropbox.com/sh/37n7drt7q4sjrzn/AAAa1HFeeRAE5M3YWG9Ac2Bta>.

After this, you can compute the trans-ethnic genetic correlation using the scores and your summary statistics. The summary statistics file must contain:
- SNP information in the form an 'rsid' or 'SNP' column with ID names or 'chr' and 'pos' columns with the chromosome number (no "chr") and base position.
- Allele information in the form of 'a1' and 'a2' or 'A1' and 'A2' columns
- 'N' column with per-SNP or whole study sample size
- 'beta' and 'SE', 'OR' and 'p-value', or 'Z' columns represing allele significance. It doesn't matter whether the effect size is signed with respect to A1 as the effect allele or A2 as the effect allele as long as it is the same in the file for both populations (eg. signed wrt A1 as effect in pop1 and A2 as effect in pop2).
- (Optional) an allele frquency 'AF' column with the frequency of A2. This will not be used except to align the effect direction of A/T or G/C SNPs which are otherwise discarded.

The command to compute the heritability and genetic correlation is then:

`popcorn fit -v 1 --cfile scores.txt --sfile1 /path/to/POP1_sumstats.txt --sfile2 /path/to/POP2_sumstats.txt correlation.txt`  

For a full list of arguments and documentations type:  

`popcorn compute -h`  
`popcorn fit -h`  

# Output:  
Popcorn reports the common-SNP observed scale (Val (Obs)) heritability (h^2) for both populations (h1^2 and h2^2),
the genetic effect or genetic impact correlation (pge or pgi), and the standard error of
these estimates (SE). If a transormation to the underlying liability scale is requested,
it also outputs the value on the liability scale (Val (Lia)). Popcorn also computes Z-scores
and p-values corresponding to these Z-scores. In the case of heritability, the p-value
reported is for a test that the heritability is *greater than 0*: P(h^2 > 0.0). In the
case of genetic correlation the p-value reported is for a test that the genetic correlation
is *less than 1.0*: P(pg < 1.0).

# Analysis considerations and notable options
As of version 1.0, the hidden `--use_regression` option has been replaced with a visible `--use_mle` option, which defaults to `False`. After extensive post-publication testing and discussion with users, we have determined that the ML estimator can be too unstable and have replaced it with a regression based fit method. This is now the default. To try with the MLE instead pass the `--use_mle` flag. We believe that regression is the appropriate method for most phenotpes.

If you have a binary phenotype and only know the number of cases and controls, you should supply `N = N_cases + N_controls` in the summary statistics file. The convestion to the liability scale (ie options `--K1`, `--K2`, `--P1`, `--P2`)will account for the imbalance between cases and controls in the study. 

Preprocessing of data summary statistics is not supported by Popcorn and must be done prior to analysis, however basic input checking is performed. Popcorn can filter on MAF and will attempt to align any mismatched alleles between the summary statistics and score files. Any other input filtering should be done prior to analysis. One example of this is the MHC regrion on Chromosome 6. If you do not filter the MHC before computing scores you might need to increase the number of SNPs in memory with the flag `--SNPs_to_store` from the default of 10000 to 20000 or more, depending on the density of your reference panel.

Another consideration is whether to compute the genetic effect correlation or the genetic impact correlation. We generally believe genetic effect is a more realistic model, but most within population analyses use a quantity close to genetic impact. To use genetic effect, pass `--gen_effect` when computing scores and fitting the model. Note that in our analysis we used a conservative MAF threshold of 0.05. We suggesting using at least 0.01. The larger the cutoff, the closer the genetic effect and genetic impact results will be.

In our analysis, we always tested the hypothesis that the genetic correlation was < 1.0 and the p-value reported is for that hypothesis. In other situations you might want to test that the genetic correlation is > 0.0.

Finally, please note that while this method can compute heritability and single population genetic correlation, the method it uses is not the same as LD score regression. The heritability and single popualtion genetic correlation estimates may differ, sometimes substantially, from the results LDscore reports, but generally the 95% confidence intervals for the parameter estimates should overlap. If you get wildly different results with the two programs, please contact us. The regression based model should give closer results to LDSR than the MLE model.

# Dependences:  
Popcorn was developed using the following external python libraries.
If you experience problems running Popcorn, try updating your libraries,
epecially less stable ones like Pandas, Pysnptools and Statsmodels,
to the versions listed here.  
numpy 1.14.2  
scipy 1.0.1  
pandas 0.22.0  
pysnptools 0.3.9  
bottleneck 1.0.0  
statsmodels 0.8.0  
(to use --plot_likelihood) matplotlib 1.5.1

# Test files:   
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
`python -m Popcorn fit -v 1 --use_mle --cfile Popcorn/test/EUR_EAS_test_ge.cscore --gen_effect --sfile1 Popcorn/test/EUR_test.txt --sfile2 Popcorn/test/EAS_test.txt Popcorn/test/EAS_EUR_test_ge_result.txt`   
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
h1^2   0.172133 NaN
h2^2   0.086708 NaN
pg     0.302372 NaN
Performing jackknife estimation of the standard error using 200 blocks of size 219 .
Jackknife iter: 200
      Val (obs)        SE          Z     P (Z)
h1^2   0.172133  0.012870  13.374388      0
h2^2   0.086708  0.008555  10.135759      0
pge    0.302372  0.053021  13.157553      0
Analysis finished at DATE
Total time elapsed: T
~~~

For more details of the method please see: http://biorxiv.org/content/early/2016/02/23/036657 or http://www.cell.com/ajhg/abstract/S0002-9297(16)30135-5
