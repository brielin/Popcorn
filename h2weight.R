# R code to generate *.h2weight files
# The *.h2weight files are used in Popcorn compute mode with --h2weight option.

# Precomputed *.h2weight files are available here
# (so you don't need to run this R code)
# http://103.253.147.127/popcorn/h2weight.20170802.zip
# http://www.fumihiko.takeuchi.name/popcorn/h2weight.20170802.zip

# The purpose of --h2weight option is to incorporate the dependence of
# per-SNP heritability on allele frequency and LD-related functional annotations.
# The *.h2weight files specify the prior weight.
#
# Allele frequency weighting is alpha = -0.25 according to
# [Speed et al. (2017) doi:10.1038/ng.3865].
# LD-related weighting is taken from
# [Gazal et al. (2017) doi:10.1038/ng.3954 Fig. 3c, Table S8a].

# The following files are required to run this R code:
# Allele frequency data of 1000G population in PLINK .frq format.
# LD-related data files baselineLD.*.annot distributed
#  https://data.broadinstitute.org/alkesgroup/LDSCORE/

for (p in c("EAS", "EUR", "SAS")) { # population
  print(p);
  
  # load input data
  c = 1 # chromosome
  data = read.table(
    paste0("~/human/1000G/phase3_shapeit2/",p,".chr",c,".maf001.snv.frq"),
    header=T, stringsAsFactors=F)
  ld = read.table(
    paste0("~/human/ldsc/1000G_Phase3_baselineLD_ldscores/baselineLD.",c,".annot_c2_73-78"),
    header=T, stringsAsFactors=F)
  ldnames = names(ld)[-1]
  print(table(data$SNP %in% ld$SNP))
  data[, ldnames] = ld[match(data$SNP, ld$SNP), ldnames]
  for (c in (2:22)) {
    foo = read.table(
      paste0("~/human/1000G/phase3_shapeit2/",p,".chr",c,".maf001.snv.frq"),
      header=T, stringsAsFactors=F)
    ld = read.table(
      paste0("~/human/ldsc/1000G_Phase3_baselineLD_ldscores/baselineLD.",c,".annot_c2_73-78"),
      header=T, stringsAsFactors=F)
    print(table(foo$SNP %in% ld$SNP))
    foo[, ldnames] = ld[match(foo$SNP, ld$SNP), ldnames]
    data = rbind(data, foo)
  }

  # compute h2weight
  taulist = c(-0.24, -0.20, -0.20, -0.13, 0.11, 0.23)
  for (i in 1:6) {
    ldname = ldnames[i];
    tau = taulist[i];
    x = data[, ldname];
    data[, ldname] = exp(tau*(x - mean(x, na.rm=T))/sd(x, na.rm=T))
  }
  #
  alpha = -0.25
  data$MAFweight = (2 * data$MAF * (1 - data$MAF))^(1 + alpha)
  data$h2weight = apply(data[, c(ldnames,"MAFweight")], 1, prod)

  # output
  #print(cor(log(data[, c(ldnames,"MAFweight","h2weight")]), use="complete.obs"))
  for (c in 1:22) {
    output = data[data$CHR==c, c("SNP", "h2weight")];
    # check again if SNPs match with plink bed/bim/fam
    foo = read.table(
      paste0("~/human/1000G/phase3_shapeit2/",p,".chr",c,".maf001.snv.bim"),
      header=F, stringsAsFactors=F)
    if(!all(output$SNP==foo$V2)) { print("ERROR"); break }
    write.table(output,
                paste0("~/human/1000G/phase3_shapeit2/",p,".chr",c,".maf001.snv.h2weight"),
                quote=F, row.names=F, col.names=F, sep="\t")
  }
}
