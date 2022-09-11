# eQTL to TWAS

This repo contains code for deriving TWAS models from eQTL summary statistics. The same code can be used for other molecular features as well (e.g. methylation, protein, chromatin). The models are in [FUSION ](http://gusevlab.org/projects/fusion/) format. 

Please refer to and cite [our paper]() comparing approaches for integration of eQTL summary statistics with GWAS summary statistics.



***

### compute_weights.R

This script uses eQTL summary statistics to create TWAS models in a format consistent with FUSION. 

The essential software and data requirements vary according to the methods selected to create models. The methods currently implemented are:

- top1 - Single strongest eQTL
- [PRS-CS](https://github.com/getian107/PRScs) - Bayesian polygenic scoring method
- [SBayesR](https://cnsgenomics.com/software/gctb/#Overview) - Bayesian polygenic scoring method, implemented in GCTB
- [DBSLMM](https://biostat0903.github.io/DBSLMM/) - Bayesian polygenic scoring method
- [LDpred2](https://privefl.github.io/bigsnpr/articles/LDpred2.html) - Bayesian polygenic scoring method, implemented in bigsnpr
- [lassosum](https://github.com/tshmak/lassosum) - Lasso-regularised polygenic scoring method
- [SuSiE](https://stephenslab.github.io/susieR/index.html) - SNP fine-mapping software, with polygenic scoring applications



<u>Essential software and data:</u>

- eQTL summary statistics

  - Must have the following columns: GENE, CHR, SNP, BP, A1, A2, BETA, SE, Z, P, FREQ, N
  - BP must be in genome build GRCh37/hg19

- Reference plink files

  - For example, 1000 Genomes Phase 3
  - Split by chromosome, containing RSIDs, build GRCh37/hg19

- R libraries

  ```{r}
  install.packages(c('optparse','data.table'))
  ```

- [GCTB software](https://cnsgenomics.com/software/gctb/#Overview) and [GCTB reference data](https://zenodo.org/record/3376628#.Yx4JsXbMKUk)
- [GCTA software](https://yanglab.westlake.edu.cn/software/gcta/#Overview) 



<u>Method specific software and data::</u>

- PRS-CS

  - [PRS-CS software](https://github.com/getian107/PRScs)
  - [PRS-CS reference data](https://github.com/getian107/PRScs#getting-started)

- DBSLMM

  - [DBSLMM software](https://biostat0903.github.io/DBSLMM/)

- LDpred2

  - [LDpred2 reference data](https://figshare.com/articles/dataset/European_LD_reference_with_blocks_/19213299)

  - ```{r}
    install.packages(c('bigsnpr','dplyr')
    ```

- lassosum

  - ```{r}
    install.packages(c('lassosum','fdrtool'))
    ```

- SuSiE

  - ```{r}
    install.packages('susieR')
    ```



<u>Example of running the script for a single gene:</u>

```{bash}
Rscript eQTL_to_TWAS/compute_weights.R \
    --id ENSG00000206503 \
    --extract /Data/ldsc/w_hm3.snplist \
    --sumstats /Data/eQTLGen/eQTLGen_sumstats.txt \
    --gcta /Software/gcta_1.94 \
    --gctb /Software/gctb_203 \
    --gctb_ref /Data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr \
    --plink_ref_chr /Data/1KG/Phase3/1KGPhase3.w_hm3.chr \
    --plink_ref_keep /Data/1KG/Phase3/keep_files/EUR_samples.keep \
    --ld_blocks /Data/LDetect/EUR \
    --rscript Rscript \
    --dbslmm /Software/DBSLMM/software \
    --plink /Software/plink1.9 \
    --PRScs_path /Software/PRScs \
    --PRScs_ref_path /Software/PRScs/ldblk_1kg_eur \
    --ldpred2_ref_dir /Data/ldpred2_ref \
    --output /Data/eQTL_to_TWAS/Test
```



This will create a file called 'ENSG00000206503.RDat' in the folder '/Data/eQTL_to_TWAS/Test'. This can be parallelised across genes, or run linearly for all genes in the --sumstats file by removing the --id parameter.



***





