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
  make_option("--models", action="store", default=c('top1','prscs','sbayesr','sbayesr_robust','susie','dbslmm','lassosum','ldpred2'), type='character',
    help="List of models to use for generating weights[optional]"),
  
  make_option("--plink_ref_chr", action="store", default=NA, type='character',
    help="Path to per chromosome refrence plink files [optional]"),
  make_option("--plink_ref_keep", action="store", default=NA, type='character',
    help="Path to keep file for plink reference [optional]"),
  
  make_option("--gctb", action="store", default=NA, type='character',
    help="Path to GCTB [required]"),
  make_option("--gctb_ref", action="store", default=NA, type='character',
    help="Path to GCTB LD reference [required]"),
  
  make_option("--PRScs_path", action="store", default=NA, type='character',
    help="Path to PRScs binary [optional]"),
  make_option("--PRScs_ref_path", action="store", default=NA, type='character',
    help="Path to PRScs LD reference [optional]"),
  
  make_option("--ld_blocks", action="store", default=NA, type='character',
    help="Path LD block data [optional]"),
  
  make_option("--rscript", action="store", default=NA, type='character',
    help="Path to Rscript binary [optional]"),
  
  make_option("--dbslmm", action="store", default=NA, type='character',
    help="Path to DBSLMM software [optional]"),
  
  make_option("--plink", action="store", default=NA, type='character',
    help="Path to PLINK v1.9 binary [optional]"),
  
  make_option("--ldpred2_ref_dir", action="store", default=NA, type='character',
              help="Path to directory containing LDPred2 reference data [optional]"),

  make_option("--sdpr", action="store", default=NA, type='character',
              help="Path to sdpr binary [optional]"),

  make_option("--hsq_p", action="store", default=NA, type='numeric',
              help="SNP-h2 p-value threshold to create TWAS models [optional]"),
  
  make_option("--output", action="store", default=NA, type='character',
    help="Name of output directory [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
system(paste0('mkdir -p ',opt$output))

models<-unlist(strsplit(opt$models, ','))

if(any(models == 'susie')){
  library(susieR)
}

if(any(models == 'lassosum')){
  library(lassosum)
  orig_wd<-getwd()
}

if(any(models == 'ldpred2')){
  library(bigsnpr)
}

cat(
  '#################################################################
# compute_weights_V2.R
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
  
  ###########
  # Estimate heritability using SBayesR
  ###########

  # Format for SBayesR
  ss_gene_i_sbayesr<-ss_gene_i[,c('SNP','A1','A2','FREQ','BETA','SE','P','N'), with=F]
  names(ss_gene_i_sbayesr)<-c('SNP','A1','A2','freq','b','se','p','N')
  
  # write in SBayesR format
  dir.create(paste0(opt$output,'/', gene_i,'/SBayesR'), recursive = T)
  fwrite(ss_gene_i_sbayesr, paste0(opt$output,'/', gene_i,'/SBayesR/',gene_i,'.txt'), sep=' ', na = "NA", quote=F)
  
  # Run SBayesR
  log<-system(paste0(opt$gctb,' --sbayes R --ldm ',opt$gctb_ref,chr_i,'.ldm.sparse --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1 --gwas-summary ',opt$output,'/', gene_i,'/SBayesR/',gene_i,'.txt --chain-length 10000 --exclude-mhc --burn-in 2000 --out-freq 1000 --out ',opt$output,'/', gene_i,'/SBayesR/',gene_i,'.SBayesR'), intern=T)
  
  # Read SbayesR heritability result
  if(file.exists(paste0(opt$output,'/', gene_i,'/SBayesR/',gene_i,'.SBayesR.parRes'))){
    
    par_res_file_i<-fread(paste0(opt$output,'/', gene_i,'/SBayesR/',gene_i,'.SBayesR.parRes'))
    par_res_file_i<-par_res_file_i[par_res_file_i$V1 == 'hsq',2:3, with=F]
    par_res_file_i$P<-pnorm(-abs(par_res_file_i$Mean/par_res_file_i$SD))
    par_res_file_i$Gene<-gene_i
    par_res_file_i<-par_res_file_i[,c('Gene','Mean','SD','P'), with=F]
    names(par_res_file_i)<-c('gene','hsq','se','p')
    
  } else {
    par_res_file_i<-data.frame(gene=gene_i,
                               hsq=NA,
                               se=NA,
                               p=NA)
  }

  # If SNP-h2 p < opt$hsq_p, generate SNP-weights
  if(is.na(opt$hsq_p) | (par_res_file_i$p < opt$hsq_p & !is.na(par_res_file_i$p))){
    
    # Create hsq object
    hsq<-c(par_res_file_i$hsq, par_res_file_i$se)
    
    # Create hsq.pv object
    hsq.pv<-par_res_file_i$p
    
    if(is.na(par_res_file_i$hsq) | par_res_file_i$hsq <= 0){
      hsq_val<-0.1
    } else {
      hsq_val<-par_res_file_i$hsq
    }
    
    # Create N.tot object
    N.tot<-max(ss_gene_i$N)
    GWAS_N<-mean(ss_gene_i$N)
    
    # Create snps object
    ss_gene_i$POS<-0
    snps<-ss_gene_i[,c('CHR','SNP','POS','BP','A1','A2')]
    names(snps)<-paste0('V',1:length(names(snps)))
    
    # Create cv.performance object
    cv.performance<-data.frame(matrix(NA, nrow=2, ncol=length(models)), row.names=c('rsq','pval'))
    names(cv.performance)<-models
    
    # Create wgt.matrix object
    wgt.matrix<-matrix(NA, ncol=length(models), nrow=nrow(ss_gene_i))
    dimnames(wgt.matrix)<-list(snps$V2,models)
    
    # Make a reference to harmonise weights from each method
    ref_tmp<-ss_gene_i[, c('SNP','A1','A2'), with=F]
    
    # Subset reference plink files
    if(!is.na(opt$plink_ref_keep)){
      cat('ref_keep used to subset reference genotype data.\n')

      dir.create(paste0(opt$output,'/', gene_i,'/ref'))
      write.table(ref_tmp$SNP, paste0(opt$output,'/', gene_i,'/ref/ref.snplist'), col.names=F, row.names = F, quote=F)
      
      system(paste0(opt$plink,' --bfile ',opt$plink_ref_chr,chr_i,' --extract ',opt$output,'/', gene_i,'/ref/ref.snplist --keep ',opt$plink_ref_keep,' --make-bed --out ',opt$output,'/', gene_i,'/ref/ref_chr',chr_i))

      system(paste0('rm ',opt$output,'/', gene_i,'/ref/ref.snplist'))
      
    } else {
      
      system(paste0(opt$plink,' --bfile ',opt$plink_ref_chr,chr_i,' --extract ',opt$output,'/', gene_i,'/ref/dbslmm_ref.snplist --make-bed --out ',opt$output,'/', gene_i,'/ref/ref_chr',chr_i))

    }
    
    for(mod in unlist(strsplit(opt$models, ','))){
      
      #####
      # top1
      #####
      if(mod == 'top1'){
        wgt.matrix[,colnames(wgt.matrix) == 'top1']<-ss_gene_i$Z
      }
      
      #####
      # sbayesr
      #####
      if(mod == 'sbayesr'){
        if(file.exists(paste0(opt$output,'/', gene_i,'/SBayesR/',gene_i,'.SBayesR.parRes'))){
          # SBayesR has already been run, so just read in the SNP-weights
          sbayesr_score<-fread(paste0(opt$output,'/', gene_i,'/SBayesR/',gene_i,'.SBayesR.snpRes'))
          sbayesr_score<-sbayesr_score[,c('Name','A1','A2','A1Effect'), with=F]
          names(sbayesr_score)<-c('SNP','A1','A2','BETA')
          
          # Flip effects so allele match eQTL sumstats
          sbayesr_score<-merge(ref_tmp, sbayesr_score, by='SNP', all=T)
          sbayesr_score$BETA[which(sbayesr_score$A1.x == sbayesr_score$A2.y)]<--sbayesr_score$BETA[which(sbayesr_score$A1.x == sbayesr_score$A2.y)]
          sbayesr_score<-sbayesr_score[,c('SNP','A1.x','A2.x','BETA'), with=F]
          names(sbayesr_score)<-c('SNP','A1','A2','BETA')
          
          # Sort score file according ss_gene_i
          sbayesr_score<-sbayesr_score[match(ss_gene_i$SNP, sbayesr_score$SNP),]
          
          wgt.matrix[,colnames(wgt.matrix) == 'sbayesr']<-sbayesr_score$BETA
        }
      }
      
      #####
      # sbayesr_robust
      #####
      if(mod == 'sbayesr_robust'){
        dir.create(paste0(opt$output,'/', gene_i,'/SBayesR_robust'))
        log<-system(paste0(opt$gctb,' --sbayes R --ldm ',opt$gctb_ref,chr_i,'.ldm.sparse --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1 --gwas-summary ',opt$output,'/', gene_i,'/SBayesR/',gene_i,'.txt --robust --chain-length 10000 --exclude-mhc --burn-in 2000 --out-freq 1000 --out ',opt$output,'/', gene_i,'/SBayesR_robust/',gene_i,'.SBayesR'), intern=T)
        
        if(file.exists(paste0(opt$output,'/', gene_i,'/SBayesR_robust/',gene_i,'.SBayesR.parRes'))){
          
          # Read in the results
          sbayesr_robust_score<-fread(paste0(opt$output,'/', gene_i,'/SBayesR_robust/',gene_i,'.SBayesR.snpRes'))
          sbayesr_robust_score<-sbayesr_robust_score[,c('Name','A1','A2','A1Effect'), with=F]
          names(sbayesr_robust_score)<-c('SNP','A1','A2','BETA')
          
          # Flip effects so allele match eQTL sumstats
          sbayesr_robust_score<-merge(ref_tmp, sbayesr_robust_score, by='SNP', all=T)
          sbayesr_robust_score$BETA[which(sbayesr_robust_score$A1.x == sbayesr_robust_score$A2.y)]<--sbayesr_robust_score$BETA[which(sbayesr_robust_score$A1.x == sbayesr_robust_score$A2.y)]
          sbayesr_robust_score<-sbayesr_robust_score[,c('SNP','A1.x','A2.x','BETA'), with=F]
          names(sbayesr_robust_score)<-c('SNP','A1','A2','BETA')
          
          # Sort score file according ss_gene_i
          sbayesr_robust_score<-sbayesr_robust_score[match(ss_gene_i$SNP, sbayesr_robust_score$SNP),]
          
          wgt.matrix[,colnames(wgt.matrix) == 'sbayesr_robust']<-sbayesr_robust_score$BETA
        }
      }
      
      #####
      # dbslmm
      #####
      if(mod == 'dbslmm'){
        # Convert to GEMMA format
        ss_gene_i_dbslmm<-ss_gene_i
        ss_gene_i_dbslmm$N_MISS<-max(ss_gene_i_dbslmm$N)-ss_gene_i_dbslmm$N
        ss_gene_i_dbslmm<-ss_gene_i_dbslmm[,c('CHR','SNP','BP','N_MISS','N','A1','A2','FREQ','BETA','SE','P'),with=F]
        names(ss_gene_i_dbslmm)<-c('chr','rs','ps','n_mis','n_obs','allele1','allele0','af','beta','se','p_wald')
        
        # Match allele1 and 0 with A1 and 2 in reference (DBSLMM calls this allele discrepancy)
        ref_bim_subset<-fread(paste0(opt$output,'/', gene_i,'/ref/ref_chr',chr_i,'.bim'))
        
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
        nsnp<-nrow(GWAS)
        
        # Write out formatted sumstats
        dir.create(paste0(opt$output,'/', gene_i,'/DBSLMM'))
        fwrite(GWAS, paste0(opt$output,'/', gene_i,'/DBSLMM/',gene_i,'.DBSLMM.txt'), sep='\t', col.names=F)
        
        # Run dbslmm
        system(paste0('chmod 777 ',opt$dbslmm,'/dbslmm'))
        system(paste0(opt$rscript,' ',opt$dbslmm,'/DBSLMM.R --plink ',opt$plink,' --block ',opt$ld_blocks,'/fourier_ls-chr',chr_i,'.bed --dbslmm ',opt$dbslmm,'/dbslmm --h2 ',hsq_val,' --ref ',opt$output,'/', gene_i,'/ref/ref_chr',chr_i,' --summary ',opt$output,'/', gene_i,'/DBSLMM/',gene_i,'.DBSLMM.txt --n ',round(GWAS_N,0),' --nsnp ',nsnp,' --outPath ',opt$output,'/', gene_i,'/DBSLMM/ --thread 1'))
        
        # Read in the results
        if(file.exists(paste0(opt$output,'/', gene_i,'/DBSLMM/',gsub('\\..*','',gene_i),'.dbslmm.txt'))){
          dbslmm_score<-fread(paste0(opt$output,'/', gene_i,'/DBSLMM/',gsub('\\..*','',gene_i),'.dbslmm.txt'))
          dbslmm_score<-dbslmm_score[,c(1,2,4), with=T]
          names(dbslmm_score)<-c('SNP','A1','BETA')
          
          # Flip effects so allele match eQTL sumstats
          dbslmm_score<-merge(ref_tmp, dbslmm_score, by='SNP', all=T)
          dbslmm_score$BETA[which(dbslmm_score$A1.x != dbslmm_score$A1.y)]<--dbslmm_score$BETA[which(dbslmm_score$A1.x != dbslmm_score$A1.y)]
          dbslmm_score<-dbslmm_score[,c('SNP','A1.x','A2','BETA'), with=F]
          names(dbslmm_score)<-c('SNP','A1','A2','BETA')
          
          # Sort score file according ss_gene_i
          dbslmm_score<-dbslmm_score[match(ss_gene_i$SNP, dbslmm_score$SNP),]
          
          wgt.matrix[,colnames(wgt.matrix) == 'dbslmm']<-dbslmm_score$BETA
        }
      }

      ######
      # PRScs
      ######
      
      if(mod == 'prscs'){
        # Format for PRScs
        ss_gene_i_prscs<-ss_gene_i[,c('SNP','A1','A2','BETA','P'), with=F]
        names(ss_gene_i_prscs)<-c('SNP','A1','A2','BETA','P')
        
        # write in PRScs format
        dir.create(paste0(opt$output,'/', gene_i,'/PRScs'))
        fwrite(ss_gene_i_prscs, paste0(opt$output,'/', gene_i,'/PRScs/',gene_i,'.txt'), sep=' ', na = "NA", quote=F)
        
        system(paste0(opt$PRScs_path,' --ref_dir=',opt$PRScs_ref_path,' --bim_prefix=',opt$output,'/', gene_i,'/ref/ref_chr',chr_i,' --sst_file=',opt$output,'/', gene_i,'/PRScs/',gene_i,'.txt --n_gwas=',round(GWAS_N,0),' --out_dir=',opt$output,'/', gene_i,'/PRScs/',gene_i,' --chrom=',chr_i))
        
        # Read in the results
        prscs_score<-fread(paste0(opt$output,'/', gene_i,'/PRScs/',gene_i,'_pst_eff_a1_b0.5_phiauto_chr',chr_i,'.txt'))
        
        skip_to_next<-F
        tryCatch(prscs_score<-prscs_score[,c('V2','V4','V6'), with=F], error = function(e){skip_to_next <<- TRUE})
        
        if(skip_to_next == F){
          names(prscs_score)<-c('SNP','A1','BETA')
          
          # Flip effects so allele match eQTL sumstats
          prscs_score<-merge(ref_tmp, prscs_score, by='SNP', all=T)
          prscs_score$BETA[which(prscs_score$A1.x != prscs_score$A1.y)]<--prscs_score$BETA[which(prscs_score$A1.x != prscs_score$A1.y)]
          prscs_score<-prscs_score[,c('SNP','A1.x','A2','BETA'), with=F]
          names(prscs_score)<-c('SNP','A1','A2','BETA')
          
          # Sort score file according ss_gene_i
          prscs_score<-prscs_score[match(ss_gene_i$SNP, prscs_score$SNP),]
          
          wgt.matrix[,colnames(wgt.matrix) == 'prscs']<-prscs_score$BETA
        }
      }
        
      ######
      # SuSiE finemapping
      ######
      
      if(mod == 'susie'){

        # Read LD estimates for eQTL sumstats
        dir.create(paste0(opt$output,'/', gene_i,'/SuSiE'))
        write.table(ss_gene_i$SNP, paste0(opt$output,'/', gene_i,'/SuSiE/',gene_i,'_snps.txt'), col.names=F, row.names=F, quote=F)
        system(paste0(opt$plink,' --bfile ',opt$output,'/', gene_i,'/ref/ref_chr',chr_i,' --extract ',opt$output,'/', gene_i,'/SuSiE/',gene_i,'_snps.txt --r square --out ',opt$output,'/', gene_i,'/SuSiE/',gene_i))
        
        ld<-as.matrix(fread(paste0(opt$output,'/', gene_i,'/SuSiE/',gene_i,'.ld')))
        
        skip_to_next<-F
        tryCatch(fitted_rss <- susie_rss(ss_gene_i$BETA/ss_gene_i$SE, n=N.tot, ld, L = 10), error = function(e){skip_to_next <<- TRUE})
        
        # Scale BETAs by PIP
        if(skip_to_next == F){
          susie_score<-data.table(SNP=ss_gene_i$SNP,
                                  A1=ss_gene_i$A1,
                                  BETA=ss_gene_i$BETA*fitted_rss$pip)

        # Flip effects so allele match eQTL sumstats
        susie_score<-merge(ref_tmp, susie_score, by='SNP', all=T)
        susie_score$BETA[which(susie_score$A1.x != susie_score$A1.y)]<--susie_score$BETA[which(susie_score$A1.x != susie_score$A1.y)]
        susie_score<-susie_score[,c('SNP','A1.x','A2','BETA'), with=F]
        names(susie_score)<-c('SNP','A1','A2','BETA')
        
        # Sort score file according ss_gene_i
        susie_score<-susie_score[match(ss_gene_i$SNP, susie_score$SNP),]
        
        wgt.matrix[,colnames(wgt.matrix) == 'susie']<-susie_score$BETA
        
        }
      }
      
      #####
      # lassosum
      #####
      
      if(mod == 'lassosum'){
        # Calculate correlation between SNP and phenotype
        # Adapt the p2cor function to allow for very small p-values
        p2cor_new<-function(z, n){
          t <- qt(pnorm(abs(z), lower.tail = FALSE, log.p = TRUE), df = GWAS_N, 
                  lower.tail = FALSE, log.p = TRUE) * sign(z)
          return(t/sqrt(n - 2 + t^2))
        }
        
        cor<-p2cor_new(z=ss_gene_i$Z, n = GWAS_N)

        # Perform lassosum to shrink effects using a range of parameters
        setwd(system.file("data", package="lassosum"))
        
        out<-lassosum.pipeline(cor=cor, chr=ss_gene_i$CHR, pos=ss_gene_i$BP, 
                               A1=ss_gene_i$A1, A2=ss_gene_i$A2,
                               ref.bfile=paste0(opt$output,'/', gene_i,'/ref/ref_chr',chr_i), 
                               LDblocks = 'EUR.hg19')
        
        setwd(orig_wd)
        
        # Perform pseudovalidation to idenitfy the best p-value threshold
        setwd(system.file("data", package="lassosum"))
        skip_to_next<-F
        tryCatch(v <- pseudovalidate(out), error = function(e){skip_to_next <<- TRUE})
        setwd(orig_wd)

        if(skip_to_next == F){
          
          # Subset the validated lassosum model
          out2 <- subset(out, s=v$best.s, lambda=v$best.lambda)
          
          # Write out a score file
          lassosum_score<-data.table(SNP=ss_gene_i$SNP[out$sumstats$order],
                                     A1=out2$sumstats$A1,
                                     BETA=out2$beta[[1]][,1])
          
          # Flip effects so allele match eQTL sumstats
          lassosum_score<-merge(ref_tmp, lassosum_score, by='SNP', all=T)
          lassosum_score$BETA[which(lassosum_score$A1.x != lassosum_score$A1.y)]<--lassosum_score$BETA[which(lassosum_score$A1.x != lassosum_score$A1.y)]
          lassosum_score<-lassosum_score[,c('SNP','A1.x','A2','BETA'), with=F]
          names(lassosum_score)<-c('SNP','A1','A2','BETA')
          
          # Sort score file according ss_gene_i
          lassosum_score<-lassosum_score[match(ss_gene_i$SNP, lassosum_score$SNP),]
          wgt.matrix[,colnames(wgt.matrix) == 'lassosum']<-lassosum_score$BETA
        }
      }

      #####
      # ldpred2
      #####

      if(mod == 'ldpred2'){
        # Attach the "bigSNP" object in R session
        if(file.exists(paste0(opt$output,'/', gene_i,'/ref/ref.bk'))){
          system(paste0('rm ',opt$output,'/', gene_i,'/ref/ref.bk'))
        }
        snp_readBed(paste0(opt$output,'/', gene_i,'/ref/ref_chr',chr_i,'.bed'), backingfile = paste0(opt$output,'/', gene_i,'/ref/ref'))
        ref <- snp_attach(paste0(opt$output,'/', gene_i,'/ref/ref.rds'))
        
        G   <- ref$genotypes
        CHR <- ref$map$chromosome
        POS <- ref$map$physical.pos
        y   <- ref$fam$affection - 1
        NCORES <- 1

        # Format sumstats for as in LDPred2 tutorial
        ss_gene_i_ldpred2<-ss_gene_i
        ss_gene_i_ldpred2<-ss_gene_i_ldpred2[,c('CHR','SNP','BP','A1','A2','BETA','SE','N','P')]
        names(ss_gene_i_ldpred2)<-c('chr','rsid','pos','a1','a0','beta','beta_se','n_eff','p')
        
        # Harmonise with the reference
        map<-readRDS(paste0(opt$ldpred2_ref_dir,'/map.rds'))
        map<-map[,c('chr','pos','a0','a1','af_UKBB','ld')]

        skip_to_next<-F
        tryCatch(info_snp <- snp_match(ss_gene_i_ldpred2, map, match.min.prop = 0), error = function(e){skip_to_next <<- TRUE})
        
        # Scale BETAs by PIP
        if(skip_to_next == F){
          
          # Perform additional suggested QC for LDPred2
          # Remove SDss<0.5???SDval or SDss>0.1+SDval or SDss<0.1 or SDval<0.05
          sd_val <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
          
          sd_y_est = median(sd_val * info_snp$beta_se * sqrt(info_snp$n_eff))
          sd_ss = with(info_snp, sd_y_est / sqrt(n_eff * beta_se^2))
  
          is_bad <-sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
          ss_gene_i_ldpred2<-info_snp[!is_bad, ]
        
          # Create sparse LD matrix
          ## indices in 'sumstats'
          ind.chr <- which(ss_gene_i_ldpred2$chr == chr_i)
          ## indices in 'map'
          ind.chr2 <- ss_gene_i_ldpred2$`_NUM_ID_`[ind.chr]
          ## indices in 'corr_chr'
          ind.chr3 <- match(ind.chr2, which(map$chr == chr_i))
            
          corr0 <- readRDS(paste0(opt$ldpred2_ref_dir,'/LD_chr', chr_i, ".rds"))[ind.chr3, ind.chr3]
            
          if(file.exists(paste0(opt$output,'/', gene_i,'/ref/LD_GW_sparse.sbk'))){
            system(paste0('rm ',opt$output,'/', gene_i,'/ref/LD_GW_sparse.sbk'))
          }
          
          skip_to_next<-F
          tryCatch(corr <- as_SFBM(corr0, paste0(opt$output,'/', gene_i,'/ref/LD_GW_sparse'), compact = TRUE), error = function(e){skip_to_next <<- TRUE})
          
          if(skip_to_next == F){
            
            # Run LDPred2-auto
            multi_auto <- snp_ldpred2_auto(corr, ss_gene_i_ldpred2, h2_init = hsq_val,
                                           vec_p_init = seq_log(1e-4, 0.9, 30),
                                           ncores = NCORES)
            
            beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
            
            # Flip beta_auto effects to match the plink reference
            names(ref$map)<-c('chr','rsid','dist','pos','a0','a1')
            tmp<-data.frame(ss_gene_i_ldpred2[,c('chr','pos','a0','a1')], beta=-1)
            info_snp_2 <- snp_match(tmp, ref$map)
            for(i in 1:ncol(beta_auto)){
              beta_auto[,i]<-beta_auto[,i]*info_snp_2$beta
            }
            
            pred_auto <- big_prodMat(G, beta_auto, ind.col = info_snp_2[["_NUM_ID_"]])
            
            sc <- apply(pred_auto, 2, sd)
            keep <- abs(sc - median(sc)) < 3 * mad(sc)
            final_beta_auto <- rowMeans(beta_auto[, keep])
            
            # compute predictions for test set
            ldpred2_score <- data.table(SNP=ss_gene_i_ldpred2$rsid, A1=ss_gene_i_ldpred2$a1, A2=ss_gene_i_ldpred2$a0, BETA = final_beta_auto)
    
            # Sort score file according ss_gene_i
            ldpred2_score<-ldpred2_score[match(ss_gene_i$SNP, ldpred2_score$SNP),]
            wgt.matrix[,colnames(wgt.matrix) == 'ldpred2']<-ldpred2_score$BETA
          }
        }
      }
      
# Note. attempted to implement SDPR and quickPRS
# SDPR would not run. No error message reported.
# QuickPRS would not converge due to not containing enough SNPs for heritability estimation.
      
################################
      
    }
    # Create RDat file for FUSION
    save(cv.performance, 
         hsq,
         hsq.pv,
         N.tot,
         snps,
         wgt.matrix,
         file = paste0(opt$output,'/',gene_i,'.RDat'))
  }
  
  # Remove temporary files
  system(paste0('rm -r ',opt$output,'/', gene_i))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')

x <- data.frame()
write.table(x, file=paste0(opt$output,'/',gene_i,'.done'), col.names=FALSE)
