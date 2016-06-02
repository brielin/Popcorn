Popcorn
======

Popcorn is a program for estimaing the correlation of causal variant effect
sizes across populations in GWAS. Estimation involves two steps:

1. Compute mode: computing cross-population scores from a reference panel,
these roughly corresponde to the similarity of the LD at each SNP across the populations
2. Fit mode: fitting the heritability and transethnic genetic correlation of
a pair of summary statistics files to the scores

Example usage:  
python -m popcorn -v 1 --out scores.txt compute --bfile1 /path/to/EUR_refpanel --bfile2 /path/to/EAS_refpanel  
python -m popcorn -v 1 --out EUR_EAS_corr.txt fit --sfile1 /path/to/EUR_sumstats.txt --sfile2 /path/to/EAS_sumstats.txt  

For a full list of options type:  
python -m popcorn -h  
python -m popcorn compute -h  
python -m popcorn fit -h  

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


Do not hesitate to contact me with questions as this is early stage software!

For more details of the method please see: http://biorxiv.org/content/early/2016/02/23/036657