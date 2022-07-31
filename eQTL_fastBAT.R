#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
    help="File containing summary statistics for molecular features [required]"),
  make_option("--id", action="store", default=NA, type='character',
    help="ID of feature [optional]"),
  make_option("--extract", action="store", default=NA, type='character',
    help="File specifying SNPs to retain [optional]"),

  make_option("--plink_ref_chr", action="store", default=NA, type='character',
    help="Path to per chromosome refrence plink files [optional]"),
  make_option("--plink_ref_keep", action="store", default=NA, type='character',
    help="Path to keep file for plink reference [optional]"),

  make_option("--chr", action="store", default=NA, type='character',
              help="Specify which chromosome to use [optional]"),
  
  make_option("--gcta", action="store", default=NA, type='character',
    help="Path to gcta [required]"),
  
  make_option("--output", action="store", default=NA, type='character',
    help="Name of output directory [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
system(paste0('mkdir -p ',opt$output))

cat(
  '#################################################################
# eQTL_fastBAT.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')

if(is.na(opt$id)){
  # Read in sumstats
  ss<-fread(opt$sumstats)
  
  cat('Sumstats contain',nrow(ss),' variant-gene associations.\n')
  cat('Sumstats contain',length(unique(ss$SNP)),' unique variants.\n')
  cat('Sumstats contain',length(unique(ss$GENE)),' unique genes\n')

} else {
  # Extract SNPs for opt$id
  ss<-fread(cmd=paste0('grep -E "GENE|',opt$id,'" ', opt$sumstats))
  ss<-ss[ss$GENE == opt$id,]
  
  cat('Computing weights for ',opt$id,'.\n')
  cat('Afer extracting gene, sumstats contain',nrow(ss),'variants.\n')
}

if(!is.na(opt$chr)){
  ss<-ss[(ss$CHR == opt$chr),]
  
  cat('Afer extracting variants, sumstats contain',nrow(ss),'variant-gene associations.\n')
  cat('Afer extracting variants, sumstats contain',length(unique(ss$SNP)),'unique variants.\n')
}

# Extract SNPs in opt$extract
if(!is.na(opt$extract)){
  extract<-fread(opt$extract)
  names(extract)[1]<-'SNP'
  extract<-extract[['SNP']]
  ss<-ss[(ss$SNP %in% extract),]
  
  cat('Afer extracting variants, sumstats contain',nrow(ss),'variant-gene associations.\n')
  cat('Afer extracting variants, sumstats contain',length(unique(ss$SNP)),'unique variants.\n')
}

# Remove variants with MAF < 1%
ss<-ss[ss$FREQ >= 0.01 & ss$FREQ <= 0.99,]

# Identify unique list of genes
genes<-unique(ss$GENE)

for(gene_i in genes){
  print(which(genes == gene_i))
  
  dir.create(paste0(opt$output,'/',gene_i))
  
  # Subset sumstats to gene (if not already)
  ss_gene_i<-ss[ss$GENE == gene_i,]
  
  # Remove genes with missing values
  ss_gene_i<-ss_gene_i[complete.cases(ss_gene_i),]
  
  # Identify chromosome number
  chr_i<-ss_gene_i$CHR[1]
  
  # Sort by chromosome and bp
  ss_gene_i<-ss_gene_i[order(ss_gene_i$CHR, ss_gene_i$BP),]
    
  # Filter SNPs to those with N > 80% of max(N)
  ss_gene_i<-ss_gene_i[ss_gene_i$N >= 0.8*max(ss_gene_i$N),]
  
  ss_gene_i<-ss_gene_i[,c('SNP','P'),with=F]
  names(ss_gene_i)<-c('SNP','p')
  
  # Set p=0 to p=minimum value that can be stored by R
  # This shouldn't matter as we are trying to see whether there is signal or not
  ss_gene_i$p[ss_gene_i$p == 0]<-.Machine$double.xmin
  
  write.table(ss_gene_i, paste0(opt$output,'/',gene_i,'/ss.txt'), col.names=T, row.names=F, quote=F)
  write.table(c(gene_i,ss_gene_i$SNP,'END'), paste0(opt$output,'/',gene_i,'/set.txt'), col.names=F, row.names=F, quote=F)
  
  # Run fastbat
  system(paste0(opt$gcta,' --bfile ',opt$plink_ref_chr,chr_i,' --keep ',opt$plink_ref_keep,' --fastBAT ',opt$output,'/',gene_i,'/ss.txt --fastBAT-set-list ',opt$output,'/',gene_i,'/set.txt --out ',opt$output,'/',gene_i,'/res --thread-num 1'))
  
  # Remove temporary files
  system(paste0('mv ',opt$output,'/',gene_i,'/res.fastbat ',opt$output,'/',gene_i,'.fastbat'))
  system(paste0('rm -r ',opt$output,'/',gene_i))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
