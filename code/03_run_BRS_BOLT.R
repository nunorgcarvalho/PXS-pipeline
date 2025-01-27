# creates and submits scripts that compute genetic correlation
library(tidyverse)
library(data.table)
source('code/00_paths.R')

tbl_Cstat <- as_tibble(fread('scratch/BRS_models/BRS_Cstats.txt'))
set_model_factors <- tbl_Cstat$model_factors %>% unique()
set_model_factors <- set_model_factors[c(2,4)] # removes cov and cov_BMI
cohorts <- tbl_Cstat$training_cohort %>% unique()
dir_pheno <- paste0(dir_repo,'scratch/cohorts/')
#loc_pheno_ALL <- paste0(dir_repo,'scratch/cohorts/BRS_cohort_ALL.txt')

dir.create('scratch/gencorrs/', showWarnings = FALSE)
dir.create('scratch/LMM/', showWarnings = FALSE)

BOLT_tbl <- tibble(cohort='',BRS1='',BRS2='')[0,]
for (cohort in cohorts) {
  for (i in 1:length(set_model_factors)) {
    for (j in i:length(set_model_factors)) {
      BOLT_tbl <- BOLT_tbl %>% add_row(
          cohort = cohort,
          BRS1 = paste0('BRS-',cohort,'-',set_model_factors[i]),
          BRS2 = paste0('BRS-',cohort,'-',set_model_factors[j]),
        )
    }
  }
}
BOLT_tbl <- BOLT_tbl %>% mutate( type=ifelse(BRS1==BRS2,'LMM','REML') )

for (i in 1:nrow(BOLT_tbl)) {
  slice <- BOLT_tbl[i,]
  loc_pheno <- paste0(dir_pheno, 'BRS_cohort_',slice$cohort,'.txt')
  
  if (slice$BRS1 != slice$BRS2) {
    phenoCol_tag <- paste0('--phenoCol ',slice$BRS1, ' --phenoCol ',slice$BRS2)
    type_tag <- '--reml --remlNoRefine'
    sbatch_oe <- paste0(slice$type,'.',slice$BRS1,'.',slice$BRS2)
    loc_script <- paste0(dir_repo,'scratch/gencorrs/',sbatch_oe,'.sh')
  } else {
    phenoCol_tag <- paste0('--phenoCol ',slice$BRS1)
    type_tag <- paste0('--lmm --verboseStats --statsFile ',dir_repo,'scratch/LMM/LMM_',slice$BRS1,'.txt')
    sbatch_oe <- paste0(slice$type,'.',slice$BRS1)
    loc_script <- paste0(dir_repo,'scratch/LMM/',sbatch_oe,'.sh')
  }
  
  
  
  script <- cat('#!/bin/sh
#SBATCH -c 20
#SBATCH -t 3-23:59
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
--statsFileBgenSnps ',paste0(dir_repo,'scratch/LMM/LMM_',slice$BRS1,'.bgen.txt'),' 
',sep='', file=loc_script)
  
}

# then run inside each script directory:
# for file in *.sh; do sbatch "$file";done




