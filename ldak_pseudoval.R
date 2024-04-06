#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default=NULL, type='character',
              help="Path to eQTL summary statistics [required]"),
  make_option("--out", action="store", default=NULL, type='character',
              help="Name of output [required]"),
  make_option("--weights", action="store", default=NULL, type='character',
              help="File listing molecular weight RDat files (must have columns WGT,ID,CHR,P0,P1) [required]"),
  make_option("--weights_dir", action="store", default=NULL, type='character',
              help="Path to directory where weight files (WGT column) are stored [required]"),
  make_option("--ref_ld_chr", action="store", default=NULL, type='character',
              help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
  make_option("--chr", action="store", default=NULL, type='numeric',
              help="Chromosome number to use [optional]"),
  make_option("--ldak", action="store", default=NULL, type='character',
              help="Path to ldak executable [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(GenoUtils)
opt$outdir<-dirname(opt$out)
system(paste0('mkdir -p ',opt$outdir))

# Create tmpdir
tmp_dir<-tempdir()

# Read in sumstats
ss<-fread(opt$sumstats)

# Read in pos file
pos<-fread(opt$weights)

# Read in reference bim
ref_bim<-NULL
for(i in ifelse(is.null(opt$chr), 1:22, opt$chr)){
  ref_bim<-rbind(ref_bim, fread(paste0(opt$ref_ld_chr,i,'.bim')))
}
names(ref_bim)<-c('CHR','SNP','POS','BP','A1','A2')

# Subset to chrosome
if(!is.null(opt$chr)){
  ss<-ss[ss$CHR == opt$chr,]
  pos<-pos[pos$CHR == opt$chr,]
  ref_bim<-ref_bim[ref_bim$CHR == opt$chr,]
}

# Identify unique list of genes
genes<-unique(ss$GENE)

# Insert IUPAC codes
ref_bim$IUPAC<-snp_iupac(ref_bim$A1, ref_bim$A2)

pseudoval_all<-NULL
obsval_all<-NULL

for(gene_i in genes){
  # Subset sumstats to gene
  ss_gene_i<-ss[ss$GENE == gene_i,]

  chr_i <- ss_gene_i$CHR[1]

  # Format for LDAK
  ss_gene_i<-data.table(
    Predictor=ss_gene_i$SNP,
    A1=ss_gene_i$A1,
    A2=ss_gene_i$A2,
    n=ss_gene_i$N,
    Z=ss_gene_i$Z
  )

  fwrite(ss_gene_i, paste0(tmp_dir,'/',gene_i,'.sumstats'), sep=' ', quote=F, na='NA')

  # Read in the SNP-weights for this gene
  load(paste0(opt$weights_dir,'/', pos$WGT[pos$ID == gene_i]))
  snps<-data.table(snps)
  names(snps)<-c('CHR','SNP','POS','BP','A1','A2')

  # Convert NA weight values to 0 to avoid error in LDAK
  wgt.matrix[is.na(wgt.matrix)]<-0

  # Create score file for ldak
  score_file <- data.table(SNP = snps$SNP,
                           A1 = snps$A1,
                           A2 = snps$A2,
                           mean = NA,
                           wgt.matrix)

  # Retain largest effect for top1 model
  score_file[-which.max(score_file[,top1]^2), 'top1'] <- 0

  # Remove variants that don't match the reference
  snps$IUPAC<-snp_iupac(snps$A1, snps$A2)
  ref_bim_i<-ref_bim[ref_bim$CHR == chr_i,]
  ref_bim_i<-merge(ref_bim_i[,c('SNP','A1','A2','IUPAC'), with=F], snps[,c('SNP','A1','A2','IUPAC'), with=F], by='SNP')
  matched<-ref_bim_i$IUPAC.x == ref_bim_i$IUPAC.y
  flipped<-detect_strand_flip(targ = ref_bim_i$IUPAC.x, ref = ref_bim_i$IUPAC.y)
  ref_bim_i<-ref_bim_i[matched | flipped,]
  score_file<-score_file[score_file$SNP %in% ref_bim_i$SNP,]

  fwrite(score_file, paste0(tmp_dir,'/',gene_i,'.score'), sep=' ', quote=F, na='NA')

  # Run ldak
  system(paste0(opt$ldak, ' --calc-scores ', tmp_dir, '/',gene_i,'.ldak --bfile ', opt$ref_ld_chr, chr_i,' --scorefile ', tmp_dir,'/',gene_i,'.score --summary ', tmp_dir,'/',gene_i,'.sumstats --power 0 --final-effects ', tmp_dir,'/',gene_i,'.score'))

  # Read in pseudoval results
  pseudoval<-fread(paste0(tmp_dir, '/',gene_i,'.ldak.cors'))
  pseudoval$V1<-colnames(wgt.matrix)
  names(pseudoval)<-c('MODEL','R')
  pseudoval$ID<-gene_i

  pseudoval_all<-rbind(pseudoval_all, pseudoval)

  # Extract cv.performance if present
  if(!is.null(cv.performance)){
    obsval<-data.table(
      MODEL=colnames(cv.performance),
      R=sqrt(cv.performance['rsq',]),
      ID=gene_i)

    obsval_all<-rbind(obsval_all, obsval)
  }
}

if(!is.null(obsval_all)){
  fwrite(obsval_all, paste0(opt$out,'.obs_cor.txt'), sep=' ', na='NA', quote=F)
}
fwrite(pseudoval_all, paste0(opt$out,'.pseudo_cor.txt'), sep=' ', na='NA', quote=F)



