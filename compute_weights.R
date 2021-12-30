#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--extract", action="store", default=NA, type='character',
    help="File specifying SNPs to retain [optional]"),
  make_option("--gctb", action="store", default=NA, type='character',
    help="Path to GCTB [required]"),
  make_option("--gctb_ref", action="store", default=NA, type='character',
    help="Path to GCTB LD reference [required]"),
  make_option("--plink_ref_chr", action="store", default=NA, type='character',
    help="Path to per chromosome refrence plink files [required]"),
  make_option("--plink_ref_keep", action="store", default=NA, type='character',
    help="Path to keep file for plink reference [required]"),
  make_option("--PRScs_path", action="store", default=NA, type='character',
    help="Path to PRScs binary [required]"),
  make_option("--PRScs_ref_path", action="store", default=NA, type='character',
    help="Path to PRScs LD reference [required]"),
  make_option("--ld_blocks", action="store", default=NA, type='character',
    help="Path LD block data [required]"),
  make_option("--rscript", action="store", default=NA, type='character',
    help="Path to Rscrpt binary [required]"),
  make_option("--dbslmm", action="store", default=NA, type='character',
    help="Path to DBSLMM software [required]"),
  make_option("--plink", action="store", default=NA, type='character',
    help="Path to PLINK v1.9 binary [required]"),
  make_option("--output", action="store", default=NA, type='character',
    help="Name of output [required]"),
  make_option("--sumstats", action="store", default=NA, type='character',
    help="File containing summary statistics for molecular features [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(susieR)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
  '#################################################################
# compute_weights.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

# Read in sumstats
ss<-fread(opt$sumstats)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Sumstats contain',nrow(ss),' variant-gene associations.\n')
cat('Sumstats contain',length(unique(ss$SNP)),' unique variants.\n')
cat('Sumstats contain',length(unique(ss$GENE)),' unique genes\n')
sink()

# Extract SNPs in opt$extract
if(!is.na(opt$extract)){
  extract<-fread(opt$extract)
  names(extract)[1]<-'SNP'
  extract<-extract[['SNP']]
  ss<-ss[(ss$SNP %in% extract),]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Afer extracting variants, sumstats contain',nrow(ss),'variant-gene associations.\n')
  cat('Afer extracting variants, sumstats contain',length(unique(ss$SNP)),'unique variants.\n')
  sink()
}

# Remove variants with MAF < 1%
ss<-ss[ss$FREQ >= 0.01 & ss$FREQ <= 0.99,]

if(!is.na(opt$plink_ref_keep)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('ref_keep used to subset reference genotype data.\n')
  sink()
  
  #####
  # Create subset of ref files
  #####
  
  for(i in 1:22){
    system(paste0(opt$plink,' --bfile ',opt$plink_ref_chr,i,' --keep ',opt$plink_ref_keep,' --make-bed --out ',opt$output_dir,'dbslmm_ref_chr',i))
  }
  
  opt$plink_ref_subset<-paste0(opt$output_dir,'dbslmm_ref_chr')
} else {
  opt$plink_ref_subset<-opt$plink_ref_chr
}

# Idnetify unique list of genes
genes<-unique(ss$GENE)

# Run for each gene seperately
dir.create(paste0(opt$output_dir,'/SBayesR'))
dir.create(paste0(opt$output_dir,'/SBayesR_robust'))
dir.create(paste0(opt$output_dir,'/DBSLMM'))
dir.create(paste0(opt$output_dir,'/PRScs'))
dir.create(paste0(opt$output_dir,'/SuSiE'))
dir.create(paste0(opt$output_dir,'/RDat_files'))

for(gene_i in genes[1:100]){
  print(which(genes == gene_i))
  ss_gene_i<-ss[ss$GENE == gene_i,]
  
  chr_i<-ss_gene_i$CHR[1]
    
  # Filter SNPs to those with N > 80% of max(N)
  ss_gene_i<-ss_gene_i[ss_gene_i$N >= 0.8*max(ss_gene_i$N),]
  
  #######
  # SBayesR
  #######
  
  # Format for SBayesR
  ss_gene_i_sbayesr<-ss_gene_i[,c('SNP','A1','A2','FREQ','BETA','SE','P','N'), with=F]
  names(ss_gene_i_sbayesr)<-c('SNP','A1','A2','freq','b','se','p','N')
  
  # write in SBayesR format
  fwrite(ss_gene_i_sbayesr, paste0(opt$output_dir,'/SBayesR/',gene_i,'.txt'), sep=' ', na = "NA", quote=F)
  
  # Run SBayesR
  log<-system(paste0(opt$gctb,' --sbayes R --ldm ',opt$gctb_ref,chr_i,'.ldm.sparse --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1 --gwas-summary ',opt$output_dir,'/SBayesR/',gene_i,'.txt --chain-length 10000 --exclude-mhc --burn-in 2000 --out-freq 1000 --out ',opt$output_dir,'/SBayesR/',gene_i,'.SBayesR'), intern=T)
  
  # Read SbayesR heritability result
  if(file.exists(paste0(opt$output_dir,'/SBayesR/',gene_i,'.SBayesR.parRes'))){
    
    par_res_file_i<-fread(paste0(opt$output_dir,'/SBayesR/',gene_i,'.SBayesR.parRes'))
    par_res_file_i<-par_res_file_i[par_res_file_i$V1 == 'hsq',2:3, with=F]
    par_res_file_i$P<-pnorm(-abs(par_res_file_i$Mean/par_res_file_i$SD))
    par_res_file_i$Gene<-gene_i
    par_res_file_i<-par_res_file_i[,c('Gene','Mean','SD','P'), with=F]
    names(par_res_file_i)<-c('gene','hsq','se','p')

    # If SNP-h2 p < 0.01, generate SNP-weights
    if(par_res_file_i$p < 0.01){
      # Flip the effect of each method to match eqtl sumstats
      ref_tmp<-ss_gene_i[, c('SNP','A1','A2'), with=F]
      
      #############
      # SBayesR
      #############
      # SBayesR has already been run, so just read in the SNP-weights
      sbayesr_score<-fread(paste0(opt$output_dir,'/SBayesR/',gene_i,'.SBayesR.snpRes'))
      sbayesr_score<-sbayesr_score[,c('Name','A1','A2','A1Effect'), with=F]
      names(sbayesr_score)<-c('SNP','A1','A2','BETA')
      
      # Flip effects so allele match eQTL sumstats
      sbayesr_score<-merge(ref_tmp, sbayesr_score, by='SNP', all=T)
      sbayesr_score$BETA[which(sbayesr_score$A1.x == sbayesr_score$A2.y)]<--sbayesr_score$BETA[which(sbayesr_score$A1.x == sbayesr_score$A2.y)]
      sbayesr_score<-sbayesr_score[,c('SNP','A1.x','A2.x','BETA'), with=F]
      names(sbayesr_score)<-c('SNP','A1','A2','BETA')
      
      #############
      # SBayesR robust
      #############
      # SBayesR format sumstats are already available
      # So just run SBayesR with robust setting
      
      log<-system(paste0(opt$gctb,' --sbayes R --ldm ',opt$gctb_ref,chr_i,'.ldm.sparse --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1 --gwas-summary ',opt$output_dir,'/SBayesR/',gene_i,'.txt --robust --chain-length 10000 --exclude-mhc --burn-in 2000 --out-freq 1000 --out ',opt$output_dir,'/SBayesR_robust/',gene_i,'.SBayesR'), intern=T)
      
      # Read in the results
      sbayesr_robust_score<-fread(paste0(opt$output_dir,'/SBayesR_robust/',gene_i,'.SBayesR.snpRes'))
      sbayesr_robust_score<-sbayesr_robust_score[,c('Name','A1','A2','A1Effect'), with=F]
      names(sbayesr_robust_score)<-c('SNP','A1','A2','BETA')
      
      # Flip effects so allele match eQTL sumstats
      sbayesr_robust_score<-merge(ref_tmp, sbayesr_robust_score, by='SNP', all=T)
      sbayesr_robust_score$BETA[which(sbayesr_robust_score$A1.x == sbayesr_robust_score$A2.y)]<--sbayesr_robust_score$BETA[which(sbayesr_robust_score$A1.x == sbayesr_robust_score$A2.y)]
      sbayesr_robust_score<-sbayesr_robust_score[,c('SNP','A1.x','A2.x','BETA'), with=F]
      names(sbayesr_robust_score)<-c('SNP','A1','A2','BETA')
      
      #############
      # DBSLMM
      #############
      
      # Convert to GEMMA format
      ss_gene_i_dbslmm<-ss_gene_i
      ss_gene_i_dbslmm$N_MISS<-max(ss_gene_i_dbslmm$N)-ss_gene_i_dbslmm$N
      ss_gene_i_dbslmm<-ss_gene_i_dbslmm[,c('CHR','SNP','BP','N_MISS','N','A1','A2','FREQ','BETA','SE','P'),with=F]
      names(ss_gene_i_dbslmm)<-c('chr','rs','ps','n_mis','n_obs','allele1','allele0','af','beta','se','p_wald')
      
      # Match allele1 and 0 with A1 and 2 in reference (DBSLMM calls this allele discrepancy)
      ref_bim_subset<-fread(paste0(opt$plink_ref_subset,chr_i,'.bim'))
      
      GWAS_match<-merge(ss_gene_i_dbslmm, ref_bim_subset[,c('V2','V5','V6'),with=F], by.x=c('rs','allele1','allele0'), by.y=c('V2','V5','V6'))
      GWAS_switch<-merge(ss_gene_i_dbslmm, ref_bim_subset[,c('V2','V5','V6'),with=F], by.x=c('rs','allele1','allele0'), by.y=c('V2','V6','V5'))
      GWAS_switch$allele_tmp<-GWAS_switch$allele0
      GWAS_switch$allele0<-GWAS_switch$allele1
      GWAS_switch$allele1<-GWAS_switch$allele_tmp
      GWAS_switch$allele_tmp<-NULL
      GWAS_switch$beta<--GWAS_switch$beta
      GWAS_switch$af<-1-GWAS_switch$af
      GWAS<-rbind(GWAS_match, GWAS_switch)
      
      GWAS<-GWAS[order(GWAS$chr, GWAS$ps),]
      GWAS<-GWAS[,c('chr','rs','ps','n_mis','n_obs','allele1','allele0','af','beta','se','p_wald'),with=F]
      GWAS_N<-mean(GWAS$n_obs)
      nsnp<-nrow(GWAS)
      
      # Write out formatted sumstats
      fwrite(GWAS, paste0(opt$output_dir,'/DBSLMM/',gene_i,'.DBSLMM.txt'), sep='\t', col.names=F)
      
      # Run dbslmm
      system(paste0('chmod 777 ',opt$dbslmm,'/dbslmm'))
      system(paste0(opt$rscript,' ',opt$dbslmm,'/DBSLMM.R --plink ',opt$plink,' --block ',opt$ld_blocks,'/fourier_ls-chr',chr_i,'.bed --dbslmm ',opt$dbslmm,'/dbslmm --h2 ',par_res_file_i$hsq,' --ref ',opt$plink_ref_subset,chr_i,' --summary ',opt$output_dir,'/DBSLMM/',gene_i,'.DBSLMM.txt --n ',round(GWAS_N,0),' --nsnp ',nsnp,' --outPath ',opt$output_dir,'/DBSLMM/ --thread 1'))
      
      # Read in the results
      if(file.exists(paste0(opt$output_dir,'/DBSLMM/',gsub('\\..*','',gene_i),'.dbslmm.txt'))){
        dbslmm_score<-fread(paste0(opt$output_dir,'/DBSLMM/',gsub('\\..*','',gene_i),'.dbslmm.txt'))
        dbslmm_score<-dbslmm_score[,c(1,2,4), with=T]
        names(dbslmm_score)<-c('SNP','A1','BETA')
        
        # Flip effects so allele match eQTL sumstats
        dbslmm_score<-merge(ref_tmp, dbslmm_score, by='SNP', all=T)
        dbslmm_score$BETA[which(dbslmm_score$A1.x != dbslmm_score$A1.y)]<--dbslmm_score$BETA[which(dbslmm_score$A1.x != dbslmm_score$A1.y)]
        dbslmm_score<-dbslmm_score[,c('SNP','A1.x','A2','BETA'), with=F]
        names(dbslmm_score)<-c('SNP','A1','A2','BETA')
      } else {
        dbslmm_score<-data.table(SNP=ss_gene_i$SNP,
                                A1=ss_gene_i$A1,
                                BETA=NA)
      }
      
      ##############
      # PRScs
      ##############
      
      # Format for PRScs
      ss_gene_i_prscs<-ss_gene_i[,c('SNP','A1','A2','BETA','P'), with=F]
      names(ss_gene_i_prscs)<-c('SNP','A1','A2','BETA','P')
      
      # write in PRScs format
      fwrite(ss_gene_i_prscs, paste0(opt$output_dir,'/PRScs/',gene_i,'.txt'), sep=' ', na = "NA", quote=F)
      
      system(paste0(opt$PRScs_path,' --ref_dir=',opt$PRScs_ref_path,' --bim_prefix=',opt$plink_ref_subset,chr_i,' --sst_file=',opt$output_dir,'/PRScs/',gene_i,'.txt --n_gwas=',round(GWAS_N,0),' --out_dir=',opt$output_dir,'/PRScs/',gene_i,' --chrom=',chr_i))
      
      # Read in the results
      prscs_score<-fread(paste0(opt$output_dir,'/PRScs/',gene_i,'_pst_eff_a1_b0.5_phiauto_chr',chr_i,'.txt'))
      prscs_score<-prscs_score[,c('V2','V4','V6'), with=F]
      names(prscs_score)<-c('SNP','A1','BETA')
      
      # Flip effects so allele match eQTL sumstats
      prscs_score<-merge(ref_tmp, prscs_score, by='SNP', all=T)
      prscs_score$BETA[which(prscs_score$A1.x != prscs_score$A1.y)]<--prscs_score$BETA[which(prscs_score$A1.x != prscs_score$A1.y)]
      prscs_score<-prscs_score[,c('SNP','A1.x','A2','BETA'), with=F]
      names(prscs_score)<-c('SNP','A1','A2','BETA')
      
      ################
      # SuSiE finemapping
      ################
      
      # Read LD estimates for eQTL sumstats
      write.table(ss_gene_i$SNP, paste0(opt$output_dir,'/SuSiE/',gene_i,'_snps.txt'), col.names=F, row.names=F, quote=F)
      system(paste0(opt$plink,' --bfile ',opt$plink_ref_subset,chr_i,' --extract ',opt$output_dir,'/SuSiE/',gene_i,'_snps.txt --r square --out ',opt$output_dir,'/SuSiE/',gene_i))
      
      ld<-as.matrix(fread(paste0(opt$output_dir,'/SuSiE/',gene_i,'.ld')))
      
      tryCatch(fitted_rss <- susie_rss(ss_gene_i$BETA/ss_gene_i$SE, ld, L = 10), error = function(e){skip_to_next <<- TRUE})
      
      if(skip_to_next == TRUE){
        susie_score<-data.table(SNP=ss_gene_i$SNP,
                                A1=ss_gene_i$A1,
                                BETA=NA)
      } else {
        susie_score<-data.table(SNP=ss_gene_i$SNP,
                                A1=ss_gene_i$A1,
                                BETA=fitted_rss$pip)
      }
      
      # Flip effects so allele match eQTL sumstats
      susie_score<-merge(ref_tmp, susie_score, by='SNP', all=T)
      susie_score$BETA[which(susie_score$A1.x != susie_score$A1.y)]<--susie_score$BETA[which(susie_score$A1.x != susie_score$A1.y)]
      susie_score<-susie_score[,c('SNP','A1.x','A2','BETA'), with=F]
      names(susie_score)<-c('SNP','A1','A2','BETA')
      
      # Create RDat file for FUSION
      cv.performance<-as.matrix(data.frame(sbayesr=c(NA,NA),
                                           sbayesr_robust=c(NA,NA),
                                           dbslmm=c(NA,NA),
                                           prscs=c(NA,NA),
                                           susie=c(NA,NA),
                                           top1=c(NA,NA), 
                                           row.names=c('rsq','pval')))
      
      # Sort eQTL data into the same order as other score files
      ss_gene_i<-ss_gene_i[match(sbayesr_score$SNP, ss_gene_i$SNP),]
      
      hsq<-c(par_res_file_i$hsq, par_res_file_i$se)
      hsq.pv<-par_res_file_i$p
      N.tot<-max(ss_gene_i$N)
      ss_gene_i$POS<-0
      snps<-ss_gene_i[,c('CHR','SNP','POS','BP','A1','A2')]
      names(snps)<-paste0('V',1:length(names(snps)))
      
      wgt.matrix<-as.matrix(data.frame(sbayesr=sbayesr_score$BETA,
                                       sbayesr_robust=sbayesr_robust_score$BETA,
                                       dbslmm=dbslmm_score$BETA,
                                       prscs=prscs_score$BETA,
                                       susie=susie_score$BETA,
                                       top1=ss_gene_i$BETA,
                                       row.names=snps$V2))
      
      save(cv.performance, 
           hsq,
           hsq.pv,
           N.tot,
           snps,
           wgt.matrix,
           file = paste0(opt$output_dir,'/RDat_files/',gene_i,'.RDat'))
    }
  }
}

# Remove temporary files
system(paste0('rm -r ',opt$output_dir,'/PRScs'))
system(paste0('rm -r ',opt$output_dir,'/DBSLMM'))
system(paste0('rm -r ',opt$output_dir,'/SBayesR'))
system(paste0('rm -r ',opt$output_dir,'/SBayesR_robust'))
system(paste0('rm -r ',opt$output_dir,'/SuSiE'))
system(paste0('rm ',opt$output_dir,'/dbslmm_ref*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
