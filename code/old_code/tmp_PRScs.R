library(tidyverse)
library(data.table)
source('code/00_paths.R')

dir_PRScs <- '/n/scratch/users/n/nur479/PRScs_stuff/'
dir_jobs <- paste0(dir_PRScs, 'jobs/')
dir_venv <- '~/group_nuno/UKBB_EUR_LDmatrix/venv_LD/'

pheno <- as_tibble(fread(paste0(dir_scratch,'cohorts/BRS_cohort_ALL.txt')))
IIDs_test <- pheno %>% filter(sample_group=='B') %>% select(IID)

loc_out <- paste0(dir_PRScs, 'test_IIDs.list', header=FALSE)




# PRS-CS ####
PRScs.py <- paste0(dir_PRScs,'PRScs.py')
chrs <- 1:22
for (chr in chrs) {
  loc_script <- paste0(dir_jobs, aes)
  
  sbatch_oe <- paste0('chr',chr,'_PRScs_script')
  file_sh <- paste0(sbatch_oe, '.sh')
  loc_script <- paste0(dir_jobs,file_sh)
  
  
  
  
  sh_list <- c(sh_list, file_sh)
  
  script <- cat('#!/bin/sh
#SBATCH -c ',CORES,'
#SBATCH -t 0-',hours,':15
#SBATCH -p short
#SBATCH --mem=',RAM_NEEDED,'G
#SBATCH -o ',sbatch_oe,'.out
#SBATCH -e ',sbatch_oe,'.err

cd ',dir_jobs,'

source ',dir_venv,'bin/activate

python ',PRScs.py,' \\
--ref_dir=',/n/scratch/users/n/nur479/PRScs_stuff/ldblk_ukbb_eur,' \\
--chr ',chr,' \\
--MAX_RAM_GB ',RAM_NEEDED,' \\
--buffer ',buffer,' \\
--MAX_WORKERS ',CORES,' \\
--dir_SNPlists ',dir_SNPlists,' \\
--dir_AWS ',dir_AWS,' \\
--dir_out ',dir_out,'
', sep='', file=loc_script)
}

# python PRScs/PRScs.py \
# --ref_dir=/n/scratch/users/n/nur479/PRScs_stuff/ldblk_ukbb_eur \
# --bim_prefix=/n/scratch/users/n/nur479/PRScs_stuff/UKBB_EUR_test \
# --sst_file=