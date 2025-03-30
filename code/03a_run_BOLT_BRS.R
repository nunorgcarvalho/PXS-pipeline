# creates and submits scripts that compute genetic correlation for BRS models
library(tidyverse)
library(data.table)
source('code/00_paths.R')


# manual table I made w/ analyses I want to run
BOLT_tbl <- as_tibble(fread('input_data/BOLT_tbl.csv')) %>%
  mutate(BRS1 = paste0('BRS-',training_cohort,'-',factors1),
         BRS2 = paste0('BRS-',training_cohort,'-',factors2),
         loc_script=NA, file=NA)
BOLT_tbl$BRS2[BOLT_tbl$factors2==''] <- NA

dir_pheno <- paste0(dir_repo,'scratch/cohorts/')

dir.create('scratch/gencorrs/', showWarnings = FALSE)
dir.create('scratch/LMM/', showWarnings = FALSE)


# makes copy of cohort_ALL but for training only ####
loc_out <- paste0(dir_pheno, 'BRS_cohort_ALL_train.txt')
if (!file.exists(loc_out)) {
  cohort_ALL <- as_tibble(fread(paste0(dir_pheno, 'BRS_cohort_ALL.txt')))
  cohort_training <- cohort_ALL %>% filter(sample_group == 'A')
  fwrite(cohort_training, loc_out, sep='\t',na="NA", quote=FALSE)
}


###
# use this to filter BOLT_tbl to the actual scripts you need to currently run
BOLT_tbl <- BOLT_tbl %>% filter(h2_cohort == 'ALL_train', type=='LMM')
###

# loops through each job and writes script for it
for (i in 1:nrow(BOLT_tbl)) {
  slice <- BOLT_tbl[i,]
  loc_pheno <- paste0(dir_pheno, 'BRS_cohort_',slice$h2_cohort,'.txt')
  
  if (slice$type == 'REML') {
    phenoCol_tag <- paste0('--phenoCol ',slice$BRS1, ' --phenoCol ',slice$BRS2)
    type_tag <- '--reml --remlNoRefine'
    sbatch_oe <- paste0(slice$type,'.',slice$h2_cohort,'.',slice$BRS1,'.',slice$BRS2)
    filename <- paste0(sbatch_oe,'.sh')
    loc_script <- paste0(dir_scratch,'gencorrs/',filename)
  } else if (slice$type == 'LMM') {
    phenoCol_tag <- paste0('--phenoCol ',slice$BRS1)
    type_tag <- paste0('--lmm --verboseStats --statsFile ',dir_repo,'scratch/LMM/LMM.',slice$h2_cohort,'.',slice$BRS1,'.txt')
    sbatch_oe <- paste0(slice$type,'.',slice$h2_cohort,'.',slice$BRS1)
    filename <- paste0(sbatch_oe,'.sh')
    loc_script <- paste0(dir_scratch,'LMM/',filename)
  }
  
  loc_statsFileBgenSnps <- paste0(dir_scratch,'LMM/LMM.',slice$h2_cohort,'.',slice$BRS1,'.bgen.txt')
  
  BOLT_tbl$file[i] <- filename
  BOLT_tbl$loc_script[i] <- loc_script
  
  
  script <- cat('#!/bin/sh
#SBATCH -c 20
#SBATCH -t 1-23:59
#SBATCH -p medium
#SBATCH --mem=100G
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
# or, run the following in the respective directory for gencorrs vs LMM

cat('for file in',paste0(BOLT_tbl$file[BOLT_tbl$type=='LMM'], collapse=' '),'; do sbatch "$file";done')
cat('for file in',paste0(BOLT_tbl$file[BOLT_tbl$type=='REML'], collapse=' '),'; do sbatch "$file";done')




