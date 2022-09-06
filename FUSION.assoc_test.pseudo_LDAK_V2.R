suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics (must have SNP and Z column headers) [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--extract_weight", action="store", default=NA, type='character',
              help="Specify ID for given SNP weight [optional]"),
  make_option("--weights", action="store", default=NA, type='character',
              help="File listing molecular weight RDat files (must have columns WGT,ID,CHR,P0,P1) [required]"),
  make_option("--weights_dir", action="store", default=NA, type='character',
              help="Path to directory where weight files (WGT column) are stored [required]"),
  make_option("--ref_ld_chr", action="store", default=NA, type='character',
              help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
  make_option("--plink1", action="store", default='plink', type='character',
              help="Path PLINK v1.9 software binary [required]"),  
  make_option("--ldak", action="store", default=NA, type='character',
              help="Path to ldak executable [required]"),
  make_option("--ldak_map", action="store", default=NA, type='character',
              help="Path to ldak map [required]"),
  make_option("--chr", action="store", default=NA, type='character',
              help="Chromosome to analyze, currently only single chromosome analyses are performed [required]")
)

library(data.table)

opt = parse_args(OptionParser(option_list=option_list))

dir.create(paste0(opt$out,'_tmp'), recursive = T)
opt$output_dir<-paste0(opt$out,'_tmp/')

# Load in summary stats
library(data.table)
sumstat = fread(opt$sumstats)

# Load in list of weights
wgtlist = fread(opt$weights)

# Subset SNP-weights in opt$extract_weight
if(!is.na(opt$extract_weight)){
  if(any(names(sumstat) == 'GENE')){
      sumstat<-sumstat[sumstat$GENE == opt$extract_weight,]
  }
  wgtlist<-wgtlist[wgtlist$ID == opt$extract_weight,]
  opt$chr = wgtlist$CHR[1]
} else {
  if(any(names(sumstat) == 'GENE')){
    setkey(sumstat, GENE)
  }
}

# Subset SNP-weights in opt$chr
if(!is.na(opt$chr)){
  sumstat<-sumstat[sumstat$CHR == opt$chr,]
  wgtlist<-wgtlist[wgtlist$CHR == opt$chr,]
}

N = nrow(wgtlist)

if(N == 0){
  cat("No models for specified chromosome\n")
  q()
}

####
# Format reference for LDAK
####

# Must be hg19 with CHR:BP IDs
system(paste0('cp ', opt$ref_ld_chr,opt$chr,'.bed ',opt$output_dir,'/ref.bed'))
system(paste0('cp ', opt$ref_ld_chr,opt$chr,'.fam ',opt$output_dir,'/ref.fam'))
system(paste0('cp ', opt$ref_ld_chr,opt$chr,'.bim ',opt$output_dir,'/ref.bim'))

# Insert genetic distances
system(paste0(opt$plink1,' --bfile ',opt$output_dir,'/ref --cm-map ',opt$ldak_map,'/genetic_map_chr@_combined_b37.txt --make-bed --out ',opt$output_dir,'/map'))
system(paste0("cat ", opt$output_dir,"/map.bim | awk '{print $2, $3}' > ", opt$output_dir,"/map.all"))
system(paste0("awk '(NR==FNR){arr[$1]=$2;next}{print $1, $2, arr[$2], $4, $5, $6}' ",opt$output_dir,"/map.all ",opt$output_dir,"/ref.bim > ",opt$output_dir,"/tmp.bim; mv ",opt$output_dir,"/tmp.bim ",opt$output_dir,"/ref.bim"))
system(paste0('rm ',opt$output_dir,'/map*'))

# For pseudovalidation, split reference into two as scores are already created (there is reference sample overlap)
system(paste0("awk < ",opt$output_dir,"/ref.fam '(NR%2==1){print $0 > \"",opt$output_dir,"/keepa\"}(NR%2==0){print $0 > \"",opt$output_dir,"/keepb\"}'"))

bim<-fread(paste0(opt$output_dir,'/ref.bim'))
bim$IUPAC[bim$V5 == 'A' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='A']<-'W'
bim$IUPAC[bim$V5 == 'C' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='C']<-'S'
bim$IUPAC[bim$V5 == 'A' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='A']<-'R'
bim$IUPAC[bim$V5 == 'C' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='C']<-'Y'
bim$IUPAC[bim$V5 == 'G' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='G']<-'K'
bim$IUPAC[bim$V5 == 'A' & bim$V6 =='C' | bim$V5 == 'C' & bim$V6 =='A']<-'M'

## For each wgt file:
pseudo_all<-NULL
for ( w in 1:nrow(wgtlist)){
  print(w)
  if(!is.na(opt$extract_weight)){
    sumstat_w<-sumstat
  } else {
    if(any(names(sumstat) == 'GENE')){
      sumstat_w<-sumstat[.(wgtlist$ID[w])]
    } else {
      sumstat_w<-sumstat
    }
  }
  
  # Format the sumstats for LDAK
  sumstat_w$IUPAC[sumstat_w$A1 == 'A' & sumstat_w$A2 =='T' | sumstat_w$A1 == 'T' & sumstat_w$A2 =='A']<-'W'
  sumstat_w$IUPAC[sumstat_w$A1 == 'C' & sumstat_w$A2 =='G' | sumstat_w$A1 == 'G' & sumstat_w$A2 =='C']<-'S'
  sumstat_w$IUPAC[sumstat_w$A1 == 'A' & sumstat_w$A2 =='G' | sumstat_w$A1 == 'G' & sumstat_w$A2 =='A']<-'R'
  sumstat_w$IUPAC[sumstat_w$A1 == 'C' & sumstat_w$A2 =='T' | sumstat_w$A1 == 'T' & sumstat_w$A2 =='C']<-'Y'
  sumstat_w$IUPAC[sumstat_w$A1 == 'G' & sumstat_w$A2 =='T' | sumstat_w$A1 == 'T' & sumstat_w$A2 =='G']<-'K'
  sumstat_w$IUPAC[sumstat_w$A1 == 'A' & sumstat_w$A2 =='C' | sumstat_w$A1 == 'C' & sumstat_w$A2 =='A']<-'M'
  
  sumstat_w<-merge(sumstat_w, bim[, c('V2','IUPAC'), with=F], by.x='SNP', by.y='V2')
  sumstat_w<-sumstat_w[sumstat_w$IUPAC.x == sumstat_w$IUPAC.y | 
                         sumstat_w$IUPAC.x == 'R' & sumstat_w$IUPAC.y == 'Y' 
                       | sumstat_w$IUPAC.x == 'Y' & sumstat_w$IUPAC.y == 'R' 
                       | sumstat_w$IUPAC.x == 'K' & sumstat_w$IUPAC.y == 'M' 
                       | sumstat_w$IUPAC.x == 'M' & sumstat_w$IUPAC.y == 'K',]
  
  sumstat_w<-sumstat_w[,c('SNP','A1','A2','N','Z')]
  names(sumstat_w)<-c('Predictor','A1','A2','n','Z')
  fwrite(sumstat_w, paste0(opt$output_dir,'GWAS_sumstats_temp.txt'), sep=' ')
  write.table(sumstat_w$Predictor, paste0(opt$output_dir,'extract.txt'), col.names=F, row.names=F, quote=F)
  
  # Format the wgt.list as a score file
  wgt.file = paste(opt$weights_dir,"/",wgtlist$WGT[w],sep='')
  load(wgt.file)
  
  wgt.matrix_dat<-data.table(dimnames(wgt.matrix)[[1]], wgt.matrix)
  score_file<-merge(snps, wgt.matrix_dat, by.x='V2', by.y='V1')
  if(any(names(score_file) == 'top1')){
    score_file$top1[which(score_file$top1 != max(score_file$top1))]<-0
    score_file$top1[duplicated(score_file$top1)]<-0
  }
  
  score_file$Centre<-NA
  
  score_file$IUPAC[score_file$V5 == 'A' & score_file$V6 =='T' | score_file$V5 == 'T' & score_file$V6 =='A']<-'W'
  score_file$IUPAC[score_file$V5 == 'C' & score_file$V6 =='G' | score_file$V5 == 'G' & score_file$V6 =='C']<-'S'
  score_file$IUPAC[score_file$V5 == 'A' & score_file$V6 =='G' | score_file$V5 == 'G' & score_file$V6 =='A']<-'R'
  score_file$IUPAC[score_file$V5 == 'C' & score_file$V6 =='T' | score_file$V5 == 'T' & score_file$V6 =='C']<-'Y'
  score_file$IUPAC[score_file$V5 == 'G' & score_file$V6 =='T' | score_file$V5 == 'T' & score_file$V6 =='G']<-'K'
  score_file$IUPAC[score_file$V5 == 'A' & score_file$V6 =='C' | score_file$V5 == 'C' & score_file$V6 =='A']<-'M'
  
  score_file<-merge(score_file, bim[, c('V2','IUPAC'), with=F], by.x='V2', by.y='V2')
  score_file<-score_file[score_file$IUPAC.x == score_file$IUPAC.y | 
                         score_file$IUPAC.x == 'R' & score_file$IUPAC.y == 'Y' 
                       | score_file$IUPAC.x == 'Y' & score_file$IUPAC.y == 'R' 
                       | score_file$IUPAC.x == 'K' & score_file$IUPAC.y == 'M' 
                       | score_file$IUPAC.x == 'M' & score_file$IUPAC.y == 'K',]
  
  score_file<-score_file[,c('V2','V5','V6','Centre', colnames(wgt.matrix))]
  names(score_file)<-c('Predictor','A1','A2','Centre', colnames(wgt.matrix))
  score_file[is.na(score_file)] <- 0
  
  write.table(score_file, paste0(opt$output_dir,'score.effects'), col.names=T, row.names=F, quote=F)
  
  #####
  # Create pseudosummaries
  #####
  
  # Create pseudo summaries
  system(paste0(opt$ldak,' --pseudo-summaries ',opt$output_dir,'GWAS_sumstats_temp.pseudo --extract ',opt$output_dir,'extract.txt --bfile ',opt$output_dir,'/ref --summary ',opt$output_dir,'GWAS_sumstats_temp.txt --training-proportion .9 --keep ',opt$output_dir,'/keepa --allow-ambiguous YES --max-threads 1'))
  
  #####
  # Evaluate models
  #####
  
  system(paste0(opt$ldak,' --calc-scores ',opt$output_dir,'/scores --bfile ',opt$output_dir,'/ref --scorefile ',opt$output_dir,'/score.effects --summary ',opt$output_dir,'GWAS_sumstats_temp.pseudo.test.summaries --power 0 --final-effects ',opt$output_dir,'/score.effects --keep ',opt$output_dir,'/keepb --allow-ambiguous YES'))
  
  pseudo<-fread(paste0(opt$output_dir,'scores.cors'))
  pseudo$V1<-colnames(wgt.matrix)
  
  names(pseudo)<-c('MODEL','R')
  pseudo$PANEL<-wgtlist$PANEL[w]
  pseudo$WGT<-wgtlist$WGT[w]
  pseudo$ID<-wgtlist$ID[w]
  
  pseudo<-pseudo[,c('PANEL','WGT','ID','MODEL','R'), with=F]
  
  pseudo_all<-rbind(pseudo_all, pseudo)
  
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats_temp.txt'))
  system(paste0('rm ',opt$output_dir,'extract.txt'))
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats_temp.pseudo.progress'))
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats_temp.pseudo.train.summaries'))
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats_temp.pseudo.test.summaries'))
  system(paste0('rm ',opt$output_dir,'score.effects'))
  system(paste0('rm ',opt$output_dir,'scores.progress'))
  system(paste0('rm ',opt$output_dir,'scores.profile'))
  system(paste0('rm ',opt$output_dir,'scores.cors'))
  system(paste0('rm ',opt$output_dir,'scores.effects.best'))
  
}

system(paste0('rm -r ',opt$output_dir))

fwrite(pseudo_all, paste0(opt$out,'.ldak.pseudo.txt'), quote=F, sep=' ', na='NA')
