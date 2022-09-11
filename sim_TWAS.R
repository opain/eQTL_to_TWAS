#!/usr/bin/Rscript
###############################################################
# Estimate transcriptome-wide significance threshold for TWAS
###############################################################

# This script calculates p-values for all features using a random phenotype, and then saves the minimum p-value. The 5th percentile of these p-values corresponds to transcriptome-wide significance.

suppressMessages(library("optparse"))

option_list = list(
make_option("--output", action="store", default=NA, type='character',
	help="output name [required]"),
make_option("--nperm", action="store", default=100, type='numeric',
	help="number of permutations [required]"),
make_option("--ncore", action="store", default=1, type='character',
	help="number of cores for parallel computing [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Open library for parallel computing
library(foreach)
library(doMC) 

# Set up parallel environment, specifying the nimber of cores to use
registerDoMC(opt$ncore)

# Load libraries for quickly reading files and regression
library(RcppEigen)
library(data.table)

# Read in predicted expression levels for all features in the FUSION 1KG reference
GeneX_all<-fread(paste0('zcat /users/k1806347/oliverpainfel/biomarkers-brc-mh/TWAS_resource/eQTL_to_TWAS/data/1kg_pred_exp/eQTLGenall_300722.eQTL/FeaturePredictions_eQTLGen.eQTL.txt.gz'),data.table=F)

# Choose the number of permutations
Nperm<-opt$nperm

# Create progress bar
pb <- txtProgressBar(min = 0, max = Nperm, style = 3)

# Record start time to see how long the permutation procedure takes
start.time <- Sys.time()
start.time

# Each loop is a permutation which generates a random phenotype, test for an association between the random phenotype and every genes. At the end it saves the minimum p-value from each permutation (i.e. null TWAS).
res<-foreach(j=1:Nperm, .combine=rbind) %dopar% {
  GeneX_all$Pheno<-rnorm(dim(GeneX_all)[1])
	
  res_i<-NULL
  for(i in 3:(dim(GeneX_all)[2]-1)){
    mod<-summary(RcppEigen::fastLm(y=GeneX_all$Pheno, X=GeneX_all[,i]))
    res_i<-rbind(res_i,data.frame(perm=j,
                                  ID=names(GeneX_all)[i],
                                  BETA=coef(mod)[1,1],
                                  SE=coef(mod)[1,2],
                                  P=coef(mod)[1,4]))
  }
  
	setTxtProgressBar(pb, j)
	print(Sys.time())
	
	res_i
	
}

# Record when permutation procedure finished and print the time taken.
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Save the minimum p value from each permutation.
fwrite(res, paste0(opt$output,'.',Nperm,'_perm.txt'), sep=' ', na='NA', quote=F)


