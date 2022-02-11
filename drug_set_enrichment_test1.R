#***********************************
# Drug-set enrichment tests
#***********************************
library(metap)
library(dplyr)
setwd("...")  ##set working directory

###load the list of drugs for which you wish to test for enrichment; e.g. can extract ATC N05A drugs for antipsychotics
file = read.table("KEGG_N05A_antipsychotics.txt",header=T) 
scz.drugs = unique( as.character(file[,1]) ) ##assuming the drug names are listed in the 1st column

#********************
# load data 
#********************
##assume we are testing the following tissues (which are included in the filenames of permutation results); can be changed
tissues = c("TW_Brain_Anterior_cingulate_cortex_BA24",
"TW_Brain_Caudate_basal_ganglia",
"TW_Brain_Cerebellar_Hemisphere",
"TW_Brain_Cerebellum",
"TW_Brain_Cortex",
"TW_Brain_Frontal_Cortex_BA9",
"TW_Brain_Hippocampus",
"TW_Brain_Hypothalamus",
"TW_Brain_Nucleus_accumbens_basal_ganglia",
"TW_Brain_Putamen_basal_ganglia")

no_tissue = length(tissues)

##load results from the permutation test comparing drug vs disease expression profiles
setwd("/Permutation_Test/Result/SCZ")
no.drugs=3478  ##no of drugs; can be changed

pval.oneSamp.t = numeric(no_tissue)
pval.twoSamp.t = numeric(no_tissue)
pval.KS = numeric(no_tissue)
mean.z.in.matchedGroup = numeric(no_tissue)

#analysis across 10 brain regions
for (i in 1:no_tissue){
  ##load results from the permutation test comparing drug vs disease expression profiles
  CMAP = read.csv( paste("AvgRank_5methods_w_permutation_HapMap_SCZ2_",tissues[i], ".csv"
                         , sep="") )         ##the name for the 1st file is AvgRank_5methods_w_permutation_HapMap_SCZ2_TW_Brain_Anterior_cingulate_cortex_BA24.csv, etc.
  
  ##find drugs in the repositioning results that match with drugs in the indication list
  ind.match.trial = grep(   paste(scz.drugs,collapse="|"), CMAP$Drug  , ignore.case=TRUE) 
  
  #**********************************************
  # regression-based gene set analysis
  #**********************************************
  pval = CMAP$perm.p
  pval[pval==1] <- 0.999
  pval[pval==0] <- 1e-4
  
  ##transform pvalues to z-scores
  zval = qnorm(pval)
  
  #self-contained test: one sample t-test (less than zero)
  fit1 = t.test(zval[ind.match.trial], alternative = "less", mu = 0, conf.level = 0.95)
  pval.oneSamp.t[i] = fit1$p.value 
  
  # competitive test: two-sample t test
  ind.notmatch.trial = setdiff( 1:nrow(CMAP), ind.match.trial) 
  fit2 = t.test( zval[ind.match.trial], zval[ind.notmatch.trial], alternative = "less")
  pval.twoSamp.t[i]  = fit2$p.value 
  
 }

#*********************************************************************************
# meta-analyses of univariate t-tests by Fisher's method and minimum p method
#*********************************************************************************
sumlog(pval.oneSamp.t)
sumlog(pval.twoSamp.t)
minimump(pval.oneSamp.t)
minimump(pval.twoSamp.t)

