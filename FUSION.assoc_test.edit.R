suppressMessages(library('plink2R'))
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
  make_option("--force_model", action="store", default=NA, type='character',
              help="Force specific predictive model to be used, no flag (default) means select most significant cross-val. Options: blup,lasso,top1,enet"),
  make_option("--all_mod", action="store", default=T, type='logical',
              help="Set to T to run TWAS using all models"),
  make_option("--max_impute", action="store", default=0.5 , type='double',
              help="Maximum fraction of SNPs allowed to be missing per gene (will be imputed using LD). [default: %default]"),
  make_option("--min_r2pred", action="store", default=0.7 , type='double',
              help="Minimum average LD-based imputation accuracy allowed for expression weight SNP Z-scores. [default: %default]"),
  make_option("--coloc_P", action="store", default=NA, type='double',
              help="P-value below which to compute COLOC statistic [Giambartolomei et al PLoS Genet 2013]\nRequires coloc library installed and --GWASN flag. [default NA/off]"),
  make_option("--GWASN", action="store", default=NA, type='integer',
              help="Total GWAS/sumstats sample size for inference of standard GWAS effect size."),
  make_option("--PANELN", action="store", default=NA, type='character',
              help="File listing sample size for each panel for inference of standard QTL effect size, cross-referenced against 'PANEL' column in weights file"),
  make_option("--chr", action="store", default=NA, type='character',
              help="Chromosome to analyze, currently only single chromosome analyses are performed [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

allele.qc = function(a1,a2,ref1,ref2) {
        a1 = toupper(a1)
        a2 = toupper(a2)
        ref1 = toupper(ref1)
        ref2 = toupper(ref2)

	ref = ref1
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip1 = flip

	ref = ref2
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip2 = flip;

	snp = list()
	snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
	snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
	snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
	snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

	return(snp)
}

# Load in summary stats
library(data.table)
sumstat = fread(opt$sumstats)

# Load in list of weights
wgtlist = fread(opt$weights)

if ( !is.na(opt$coloc_P) ) {
  if ( is.na(opt$GWASN) || opt$GWASN < 1 ) {
    cat("ERROR : --GWASN flag required to be positive integer for COLOC analysis\n")
    q()
  }
  if ( sum(names(wgtlist) == "N") == 0 ) {
    if ( sum(names(wgtlist) == "PANEL") == 0 || is.na(opt$PANELN) ) {
      cat("ERROR : 'N' field needed in weights file or 'PANEL' column and --PANELN flag required for COLOC analysis\n")
      q()
    } else {
      paneln = read.table(opt$PANELN,as.is=T,head=T,sep='\t')
      m = match( wgtlist$PANEL , paneln$PANEL )
      wgtlist$N = paneln$N[ m ]
    }
  }
  suppressMessages(library('coloc'))
}

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
  if(any(names(sumstat) == 'CHR')){
    sumstat<-sumstat[sumstat$CHR == opt$chr,]
  }
  wgtlist<-wgtlist[wgtlist$CHR == opt$chr,]
}

N = nrow(wgtlist)

if(N == 0){
  cat("No models for specified chromosome\n")
  q()
}

# Load in reference data
genos = read_plink(paste(opt$ref_ld_chr,opt$chr,sep=''),impute="avg")

## For each wgt file:
FAIL.ctr<-0
out.tbl.all.all<-NULL
for ( w in 1:nrow(wgtlist)){
  if(!is.na(opt$extract_weight)){
    sumstat_w<-sumstat
  } else {
    if(any(names(sumstat) == 'GENE')){
      sumstat_w<-sumstat[.(wgtlist$ID[w])]
    } else {
      sumstat_w<-sumstat
    }
  }

  genos_w<-genos

  # Match summary data to input, record NA where summary data is missing
  m = match( genos_w$bim[,2] , sumstat_w$SNP )
  sum.missing = is.na(m)
  sumstat_w = sumstat_w[m,]
  sumstat_w$SNP = genos_w$bim[,2]
  sumstat_w$A1[ sum.missing ] = genos_w$bim[sum.missing,5]
  sumstat_w$A2[ sum.missing ] = genos_w$bim[sum.missing,6]

  # QC / allele-flip the input and output
  qc = allele.qc( sumstat_w$A1 , sumstat_w$A2 , genos_w$bim[,5] , genos_w$bim[,6] )

  # Flip Z-scores for mismatching alleles
  sumstat_w$Z[ qc$flip ] = -1 * sumstat_w$Z[ qc$flip ]
  sumstat_w$A1[ qc$flip ] = genos_w$bim[qc$flip,5]
  sumstat_w$A2[ qc$flip ] = genos_w$bim[qc$flip,6]

  # Remove strand ambiguous SNPs (if any)
  if ( sum(!qc$keep) > 0 ) {
  	genos_w$bim = genos_w$bim[qc$keep,]
  	genos_w$bed = genos_w$bed[,qc$keep]
  	sumstat_w = sumstat_w[qc$keep,]
  }

	# Load weights
	wgt.file = paste(opt$weights_dir,"/",wgtlist$WGT[w],sep='')
	load(wgt.file)
	# Remove NAs (these should not be here)
	wgt.matrix[is.na(wgt.matrix)] = 0

	# Match up the SNPs and weights
	snps<-as.matrix(snps)
	m = match( snps[,2] , genos_w$bim[,2] )
	m.keep = !is.na(m)
	snps = snps[m.keep,]
	if(is.matrix(snps) == F){
	  snps<-t(as.matrix(snps))
	}

	wgt.matrix = wgt.matrix[m.keep,,drop=F]
	cur.genos_w = scale(genos_w$bed[,m[m.keep]])
	cur.bim = genos_w$bim[m[m.keep],]

	# Flip WEIGHTS for mismatching alleles
	qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
	wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]

	# Match up the SNPs and the summary stats
	m = match(cur.bim[,2] , sumstat_w$SNP)
	cur.Z = sumstat_w$Z[m]

	# which rows have rsq
	row.rsq = grep( "rsq" , rownames(cv.performance) )
	# which rows have p-values
	row.pval = grep( "pval" , rownames(cv.performance) )

	if(opt$all_mod == T){
	  models<-1:ncol(wgt.matrix)
	}
	if(!is.na(opt$force_model)){
	  models<-which( colnames(wgt.matrix) == opt$force_model )
	}
	if(is.na(opt$force_model) & opt$all_mod == F){
	  models = which.min(apply(cv.performance[row.pval,,drop=F],2,min,na.rm=T))
	}

	if ( length(models) == 0 ) {
		cat( "WARNING : " , unlist(wgtlist[w,]) , " did not have a predictive model ... skipping entirely\n" )
		FAIL.ctr = FAIL.ctr + 1
		next
	}

	out.tbl.all<-NULL
	for(mod.best in models){
	  cur.FAIL = FALSE

	  # Store top1 for coloc

	  if(is.na(opt$coloc_P)){
	    out.tbl = data.frame( "PANEL" = rep(NA,1) , "FILE" = character(1) , "ID" = character(1) , "CHR" = numeric(1) , "P0" = character(1) , "P1" = character(1) ,"HSQ" = numeric(1) , "NSNP" = numeric(1) , "NWGT" = numeric(1) , "MODEL" = character(1), "MODELCV.R2" = character(1) , "TWAS.Z" = numeric(1) , "TWAS.P" = numeric(1) , stringsAsFactors=FALSE )
	  } else {
	    out.tbl = data.frame( "PANEL" = rep(NA,1) , "FILE" = character(1) , "ID" = character(1) , "CHR" = numeric(1) , "P0" = character(1) , "P1" = character(1) ,"HSQ" = numeric(1) , "NSNP" = numeric(1) , "NWGT" = numeric(1) , "MODEL" = character(1), "MODELCV.R2" = character(1) , "TWAS.Z" = numeric(1) , "TWAS.P" = numeric(1) , "COLOC.PP0" = numeric(1) , "COLOC.PP1" = numeric(1) , "COLOC.PP2" = numeric(1) , "COLOC.PP3" = numeric(1), "COLOC.PP4" = numeric(1) , stringsAsFactors=FALSE )
	  }

	 	if ( sum(wgt.matrix[, mod.best] != 0) == 0 ) {
  		cat( "WARNING : " , unlist(wgtlist[w,]) , names(cv.performance)[ mod.best ] , "had", length(cur.Z) , "overlapping SNPs, but none with non-zero expression weights, try more SNPS or a different model\n")
  		cur.FAIL = TRUE
  	}

  	# if this is a top1 model, clear out all the other weights
  	if ( substr( (colnames(cv.performance))[ mod.best ],1,4) == "top1" ) wgt.matrix[ -which.max(wgt.matrix[,mod.best]^2)  , mod.best] = 0

  	# Compute LD matrix
  	if ( length(cur.Z) == 0 ) {
  		cat( "WARNING : " , unlist(wgtlist[w,]) , " had no overlapping SNPs\n")
  		cur.FAIL = TRUE
  		out.tbl$NSNP[1] = NA
  	} else if ( !cur.FAIL ) {
  		cur.LD = t(cur.genos_w) %*% cur.genos_w / (nrow(cur.genos_w)-1)
  		out.tbl$NSNP[1] = nrow(cur.LD)
  		cur.miss = is.na(cur.Z)
  		# Impute missing Z-scores
  		if ( sum(cur.miss) != 0 ) {
  			if ( sum(!cur.miss) == 0 ) {
  				cat( "WARNING : " , unlist(wgtlist[w,]) , "had no overlapping GWAS Z-scores, skipping this gene\n")
  				cur.FAIL = TRUE
  			} else if ( mean(cur.miss) > opt$max_impute ) {
  				cat( "WARNING : " , unlist(wgtlist[w,]) , "had" , sum(cur.miss) , "/" , length(cur.miss) , "non-overlapping GWAS Z-scores, skipping this gene.\n")
  				cur.FAIL = TRUE
  			} else {
  				cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
  				cur.impz = cur.wgt %*% cur.Z[!cur.miss]
  				cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
  				cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)

  				all.r2pred = rep(1,length(cur.Z))
  				all.r2pred[ cur.miss ] = cur.r2pred
  				if ( sum(is.na(all.r2pred)) != 0 ) {
  					cat( "WARNING : " , unlist(wgtlist[w,]) , "had missing GWAS Z-scores that could not be imputed, skipping this gene.\n" )
  					cur.FAIL = TRUE
  				} else if ( mean( all.r2pred[ wgt.matrix[,mod.best] != 0 ] ) < opt$min_r2pred ) {
  					cat( "WARNING : " , unlist(wgtlist[w,]) , "had mean GWAS Z-score imputation r2 of" , mean( all.r2pred[ wgt.matrix[,mod.best] != 0 ] ) , "at expression weight SNPs, skipping this gene.\n")
  					cur.FAIL = TRUE
  				}
  			}
  		}

  		if ( !cur.FAIL ) {
  			# Compute TWAS Z-score
  			cur.twasz = wgt.matrix[,mod.best] %*% cur.Z
  			cur.twasr2pred = wgt.matrix[,mod.best] %*% cur.LD %*% wgt.matrix[,mod.best]

  			if ( cur.twasr2pred > 0 ) {
  				cur.twas = cur.twasz / sqrt(cur.twasr2pred)
  			} else {
  				cur.FAIL=T
  				cat( "WARNING : " , unlist(wgtlist[w,]) , " had zero predictive accuracy, try a different model.\n")
  			}
  		}
  	}

  	# populate the output
  	if ( sum(names(wgtlist) == "PANEL") == 1 ) out.tbl$PANEL[1] = wgtlist$PANEL[w]
  	out.tbl$FILE[1] = wgt.file
  	out.tbl$CHR[1] = wgtlist$CHR[w]
  	out.tbl$P0[1] = wgtlist$P0[w]
  	out.tbl$P1[1] = wgtlist$P1[w]
  	out.tbl$ID[1] = wgtlist$ID[w]
  	if ( exists("hsq") ) {
  		out.tbl$HSQ[1] = hsq[1]
  	}
  	out.tbl$MODEL[1] = colnames( cv.performance )[ mod.best ]
  	out.tbl$MODELCV.R2[1] = paste(format(cv.performance[row.rsq,mod.best],digits=2,trim=T),collapse=',')

  	if ( !cur.FAIL ) {
  		out.tbl$NWGT[1] = sum( wgt.matrix[,mod.best] != 0 )
  		out.tbl$TWAS.Z[1] = cur.twas
  		out.tbl$TWAS.P[1] = 2*(pnorm( abs(out.tbl$TWAS.Z[1]) , lower.tail=F))
  	} else {
  	  out.tbl$TWAS.Z[1] = NA
  	  out.tbl$TWAS.P[1] = NA
  	}

  	if ( cur.FAIL ) FAIL.ctr = FAIL.ctr + 1
  	out.tbl.all<-rbind(out.tbl.all, out.tbl)
	}

	# perform COLOC test
	if ( !is.na(opt$coloc_P) && any(!is.na(out.tbl.all$TWAS.Z[out.tbl.all$FILE == wgt.file])) && any(out.tbl.all$TWAS.P[out.tbl.all$FILE == wgt.file] < opt$coloc_P, na.rm = T) && !is.na(wgtlist$N[w]) ) {
	  b1 = wgt.matrix[,'top1'] / sqrt(wgtlist$N[w])
	  b2 = cur.Z / sqrt(opt$GWASN)

	  vb1 = rep(1/wgtlist$N[w],length(b1))
	  vb2 = rep(1/opt$GWASN,length(b2))

	  err = suppressMessages(capture.output(clc <- coloc.abf(dataset1=list(beta=b1,varbeta=vb1,type="quant",N=wgtlist$N[w],sdY=1),dataset2=list(beta=b2,varbeta=vb2,type="quant",N=opt$GWASN,sdY=1))))

	  out.tbl.all$COLOC.PP0[out.tbl.all$FILE == wgt.file]<-round(clc$summary[2],3)
	  out.tbl.all$COLOC.PP1[out.tbl.all$FILE == wgt.file]<-round(clc$summary[3],3)
	  out.tbl.all$COLOC.PP2[out.tbl.all$FILE == wgt.file]<-round(clc$summary[4],3)
	  out.tbl.all$COLOC.PP3[out.tbl.all$FILE == wgt.file]<-round(clc$summary[5],3)
	  out.tbl.all$COLOC.PP4[out.tbl.all$FILE == wgt.file]<-round(clc$summary[6],3)
	} else {
	  out.tbl.all$COLOC.PP0[out.tbl.all$FILE == wgt.file]<-NA
	  out.tbl.all$COLOC.PP1[out.tbl.all$FILE == wgt.file]<-NA
	  out.tbl.all$COLOC.PP2[out.tbl.all$FILE == wgt.file]<-NA
	  out.tbl.all$COLOC.PP3[out.tbl.all$FILE == wgt.file]<-NA
	  out.tbl.all$COLOC.PP4[out.tbl.all$FILE == wgt.file]<-NA
	}

	out.tbl.all.all<-rbind(out.tbl.all.all, out.tbl.all)
}

cat("Analysis completed.\n")
cat("NOTE:",FAIL.ctr,"genes-model pairs were skipped\n")

# WRITE MHC TO SEPARATE FILE
mhc = as.numeric(out.tbl$CHR) == 6 & as.numeric(out.tbl$P0) > 26e6 & as.numeric(out.tbl$P1) < 34e6

out.tbl.all.all$P0 = apply( as.matrix(out.tbl.all.all$P0) , 1 , toString )
out.tbl.all.all$P1 = apply( as.matrix(out.tbl.all.all$P1) , 1 , toString )

if ( sum( mhc ) > 0 ) {
	cat("Results in the MHC are written to",paste(opt$out,".MHC",sep=''),", evaluate with caution due to complex LD structure\n")
	write.table( format( out.tbl.all.all[mhc,] , digits=3 ) , quote=F , row.names=F , sep='\t' , file=paste(opt$out,".MHC",sep='') )
}
write.table( format( out.tbl.all.all[!mhc,] , digits=3 ) , quote=F , row.names=F , sep='\t' , file=opt$out )
