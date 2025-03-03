# creates and submits scripts that compute genetic correlation for CRS and CRFs
library(tidyverse)
library(data.table)
source('code/00_paths.R')

# reads relevant terms 
fields <- as_tibble(fread(paste0(dir_scratch,'fields_tbl.txt')))
CRFs <- fields$term[fields$use_type == 'CRF']
loc_pheno <- paste0(dir_scratch,'cohorts/BRS_cohort_ALL.txt')
BRS_coeffs <- as_tibble(fread(paste0(dir_scratch,'BRS_models/BRS_coefficients.txt'))) %>%
  filter(term %in% fields$term[fields$use_type == 'behavior'],
         `beta-BRS-ALL-cov_bvr` != 0)
# makes table of scripts to run
BOLT_tbl <- tibble(
    term1 = 'BRS-ALL-cov_bvr',term2 = c(CRFs, 'CRS-ALL')
  ) %>% add_row(
    term1 = 'CRS-ALL', term2 = BRS_coeffs$term) %>%
  add_row(term1 = c('CRS-ALL', CRFs) ) %>%
  mutate(type = ifelse(is.na(term2),'LMM','REML'),
         term1_ = str_replace(term1, '\\.','_'),
         term2_ = str_replace(term2, '\\.','_'),
         loc_script=NA, file=NA)

# loops through each job, writes and saves script
for (i in 1:nrow(BOLT_tbl)) {
  slice <- BOLT_tbl[i,]
  
  if (slice$type == 'REML') {
    phenoCol_tag <- paste0('--phenoCol ',slice$term1, ' --phenoCol ',slice$term2)
    type_tag <- '--reml --remlNoRefine'
    sbatch_oe <- paste0(slice$type,'.',slice$term1_,'.',slice$term2_)
    filename <- paste0(sbatch_oe,'.sh')
    loc_script <- paste0(dir_repo,'scratch/gencorrs/',filename)
  } else if (slice$type == 'LMM') {
    phenoCol_tag <- paste0('--phenoCol ',slice$term1)
    type_tag <- paste0('--lmm --verboseStats --statsFile ',dir_repo,'scratch/LMM/LMM.',slice$term1_,'.txt')
    sbatch_oe <- paste0(slice$type,'.',slice$term1_)
    filename <- paste0(sbatch_oe,'.sh')
    loc_script <- paste0(dir_repo,'scratch/LMM/',filename)
  }
  
  loc_statsFileBgenSnps <- paste0(dir_scratch,'LMM/LMM.',slice$term1_,'.bgen.txt')
  
  BOLT_tbl$file[i] <- filename
  BOLT_tbl$loc_script[i] <- loc_script
  
  script <- cat('#!/bin/sh
#SBATCH -c 20
#SBATCH -t 1-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o ',sbatch_oe,'.out
#SBATCH -e ',sbatch_oe,'.err

~/bolt \\
--numThreads 20 \\
--bed /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \\
--bim /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \\
--fam /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \\
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \\
--remove ',dir_repo,'input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \\
--phenoFile ',loc_pheno,' \\
',phenoCol_tag,' \\
--covarFile ',loc_pheno,' \\
--covarCol sex \\
--qCovarCol age \\
--qCovarCol pc{1:40} \\
--covarMaxLevels 25 \\
',type_tag,' \\
--bgenFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb_imp_chr{1:22}_v3.bgen \\
--bgenMinMAF 1e-3 \\
--bgenMinINFO 0.3 \\
--sampleFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb22881_imp_chr1_v3_s487324.sample \\
--statsFileBgenSnps ',loc_statsFileBgenSnps,' 
',sep='', file=loc_script)
  
}
# then run inside each script directory:
# for file in *.sh; do sbatch "$file";done
# or, run the following:

cat('for file in ',paste0(BOLT_tbl$file[BOLT_tbl$type=='LMM'], collapse=' '),'; do sbatch "$file";done')
cat('for file in ',paste0(BOLT_tbl$file[BOLT_tbl$type=='REML'], collapse=' '),'; do sbatch "$file";done')
