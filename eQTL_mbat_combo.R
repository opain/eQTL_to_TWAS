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
    help="Path and prefix of results file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
system(paste0('mkdir -p ',opt$output,'_tmp'))
opt$output_dir<-paste0(opt$output,'_tmp')

cat(
  '#################################################################
# eQTL_mbat_combo.R
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

# Read in the reference data
bim<-NULL
if(is.na(opt$chr)){
  for(i in 1:22){
    bim<-rbind(bim, fread(paste0(opt$plink_ref_chr,i,'.bim')))
  }
} else {
  bim<-rbind(bim, fread(paste0(opt$plink_ref_chr,opt$chr,'.bim')))
}

# Identify unique list of genes
genes<-unique(ss$GENE)

res_all<-NULL
for(gene_i in genes){
  print(which(genes == gene_i))
  
  dir.create(paste0(opt$output_dir,'/',gene_i))
  
  # Subset sumstats to gene (if not already)
  ss_gene_i<-ss[ss$GENE == gene_i,]
  
  # Remove genes with missing values
  ss_gene_i<-ss_gene_i[complete.cases(ss_gene_i),]
  
  # Identify chromosome number
  chr_i<-ss_gene_i$CHR[1]
  
  # Identify min and max BP position so fastbat includes all eQTL results.
  # Use LD reference to ensure correct build is used.
  min_bp<-min(bim$V4[bim$V2 %in% ss_gene_i$SNP])
  max_bp<-max(bim$V4[bim$V2 %in% ss_gene_i$SNP])

  # Filter SNPs to those with N > 80% of max(N)
  ss_gene_i<-ss_gene_i[ss_gene_i$N >= 0.8*max(ss_gene_i$N),]
  
  ss_gene_i<-ss_gene_i[,c('SNP','A1','A2','FREQ','BETA','SE','P','N'),with=F]
  names(ss_gene_i)<-c('SNP','A1','A2','freq','BETA','SE','P','N')
  
  # Set p=0 to p=minimum value that can be stored by R
  # This shouldn't matter as we are trying to see whether there is signal or not
  ss_gene_i$P[ss_gene_i$P == 0]<-.Machine$double.xmin
  
  # Make a gene list file
  gene_list<-data.frame( V1=chr_i,
                         V2=min_bp,
                         V3=max_bp,
                         V4=gene_i)
  
  write.table(ss_gene_i, paste0(opt$output_dir,'/',gene_i,'/ss.txt'), col.names=T, row.names=F, quote=F)
  write.table(gene_list, paste0(opt$output_dir,'/',gene_i,'/gene_list.txt'), col.names=F, row.names=F, quote=F)
  write.table(c(gene_i,ss_gene_i$SNP,'END'), paste0(opt$output_dir,'/',gene_i,'/set.txt'), col.names=F, row.names=F, quote=F)
  
  # Run mbat-combo
  system(paste0(opt$gcta,' --bfile ',opt$plink_ref_chr,chr_i,' --keep ',opt$plink_ref_keep,' --mBAT-combo ',opt$output_dir,'/',gene_i,'/ss.txt --mBAT-print-all-p --mBAT-gene-list ',opt$output_dir,'/',gene_i,'/gene_list.txt --out ',opt$output_dir,'/',gene_i,'/res --thread-num 1'))

  res<-fread(paste0(opt$output_dir,'/',gene_i,'/res.gene.assoc.mbat'))
  
  if(nrow(res) == 0){
    # For genes with only one SNP available, insert the min P
    
    ss_gene_i<-ss_gene_i[ss_gene_i$SNP %in% bim$V2,]
    bim_gene_i<-bim[bim$V2 %in% ss_gene_i$SNP,]
    
    ss_gene_i$IUPAC[ss_gene_i$A1 == 'A' & ss_gene_i$A2 =='T' | ss_gene_i$A1 == 'T' & ss_gene_i$A2 =='A']<-'W'
    ss_gene_i$IUPAC[ss_gene_i$A1 == 'C' & ss_gene_i$A2 =='G' | ss_gene_i$A1 == 'G' & ss_gene_i$A2 =='C']<-'S'
    ss_gene_i$IUPAC[ss_gene_i$A1 == 'A' & ss_gene_i$A2 =='G' | ss_gene_i$A1 == 'G' & ss_gene_i$A2 =='A']<-'R'
    ss_gene_i$IUPAC[ss_gene_i$A1 == 'C' & ss_gene_i$A2 =='T' | ss_gene_i$A1 == 'T' & ss_gene_i$A2 =='C']<-'Y'
    ss_gene_i$IUPAC[ss_gene_i$A1 == 'G' & ss_gene_i$A2 =='T' | ss_gene_i$A1 == 'T' & ss_gene_i$A2 =='G']<-'K'
    ss_gene_i$IUPAC[ss_gene_i$A1 == 'A' & ss_gene_i$A2 =='C' | ss_gene_i$A1 == 'C' & ss_gene_i$A2 =='A']<-'M'

    bim_gene_i$IUPAC[bim_gene_i$V5 == 'A' & bim_gene_i$V6 =='T' | bim_gene_i$V5 == 'T' & bim_gene_i$V6 =='A']<-'W'
    bim_gene_i$IUPAC[bim_gene_i$V5 == 'C' & bim_gene_i$V6 =='G' | bim_gene_i$V5 == 'G' & bim_gene_i$V6 =='C']<-'S'
    bim_gene_i$IUPAC[bim_gene_i$V5 == 'A' & bim_gene_i$V6 =='G' | bim_gene_i$V5 == 'G' & bim_gene_i$V6 =='A']<-'R'
    bim_gene_i$IUPAC[bim_gene_i$V5 == 'C' & bim_gene_i$V6 =='T' | bim_gene_i$V5 == 'T' & bim_gene_i$V6 =='C']<-'Y'
    bim_gene_i$IUPAC[bim_gene_i$V5 == 'G' & bim_gene_i$V6 =='T' | bim_gene_i$V5 == 'T' & bim_gene_i$V6 =='G']<-'K'
    bim_gene_i$IUPAC[bim_gene_i$V5 == 'A' & bim_gene_i$V6 =='C' | bim_gene_i$V5 == 'C' & bim_gene_i$V6 =='A']<-'M'
    
    ss_gene_i<-merge(ss_gene_i, bim_gene_i[,c('V2','IUPAC'), with=F], by.x='SNP', by.y='V2')
    ss_gene_i<-ss_gene_i[ss_gene_i$IUPAC.x == ss_gene_i$IUPAC.y | 
                           ss_gene_i$IUPAC.x == 'R' & ss_gene_i$IUPAC.y == 'Y' | 
                           ss_gene_i$IUPAC.x == 'Y' & ss_gene_i$IUPAC.y == 'R' | 
                           ss_gene_i$IUPAC.x == 'K' & ss_gene_i$IUPAC.y == 'M' | 
                           ss_gene_i$IUPAC.x == 'M' & ss_gene_i$IUPAC.y == 'K',]
    
    res_tmp<-data.frame(matrix(NA, nrow=1, ncol=ncol(res)))
    names(res_tmp)<-names(res)
    res_tmp$Gene<-gene_i
    res_tmp$Chr<-chr_i
    res_tmp$Start<-min_bp
    res_tmp$End<-max_bp
    res_tmp$No.SNPs<-1
    res_tmp$SNP_start<-min_bp
    res_tmp$SNP_end<-max_bp
    res_tmp$TopSNP<-ss_gene_i$SNP[ss_gene_i$P == min(ss_gene_i$P)][1]
    res_tmp$TopSNP_Pvalue<-min(ss_gene_i$P)
    
    res<-res_tmp
  }

  res_all<-rbind(res_all, res)
  
  # Remove temporary files
  system(paste0('rm -r ',opt$output_dir, '/',gene_i))
  
}

fwrite(res_all, paste0(opt$output, '.res.txt'), quote=F, sep=' ', na='NA')

# Remove temporary files
system(paste0('rm -r ',opt$output_dir))

end.time <- Sys.time()
time.taken <- end.time - start.time
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
