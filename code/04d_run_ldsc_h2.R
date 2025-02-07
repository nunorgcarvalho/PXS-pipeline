# computes heritability from BOLT GWASs using LDsc

# libraries and directories ####
library(tidyverse)
library(data.table)
source('code/00_paths.R')

dir_LMM <- paste0(dir_scratch,'LMM/')
dir_ldsc_out <- paste0(dir_scratch,'ldsc/') # where files and scripts will go
dir.create(dir_ldsc_out, showWarnings = FALSE)
dir_ldsc <- '~/group_nuno/ldsc/' # ldsc software
dir_LDscore <- '~/group_nuno/LDscore/' # LD scores

LMM_files <- list.files(paste0(dir_LMM), pattern='*.bgen.txt')

for (LMM_file in LMM_files) {
  filepath <- paste0(dir_LMM, LMM_file)
  GWAS_name <- str_replace(str_sub(LMM_file,5),'.bgen.txt','')
  sbatch_oe <- paste0('ldsc.',GWAS_name)
  loc_script <- paste0(dir_ldsc_out, sbatch_oe, '.sh')
  loc_out_stem <- paste0(dir_ldsc_out,'ldsc.',GWAS_name)
  
  # extracts sample size from .out file
  loc_LMM.out <- str_replace( filepath, '.bgen.txt','.out')
  N <- read_lines(loc_LMM.out) %>%
    str_subset("Number of individuals used in analysis: Nused =") %>%
    tail(1) %>% str_extract("(?<=Nused = )\\d+") %>% as.numeric()
  
  script <- cat('#!/bin/sh
#SBATCH -c 4
#SBATCH -t 0-00:30
#SBATCH -p short
#SBATCH --mem=4G
#SBATCH -o ',sbatch_oe,'.out
#SBATCH -e ',sbatch_oe,'.err

cd ',dir_ldsc_out,'

module load miniconda3
conda env create --file ',dir_ldsc,'environment.yml
source activate ldsc

echo "Formatting summary file for ldsc"
',dir_ldsc,'munge_sumstats.py \\
--sumstats ',filepath,' \\
--N ',N,' \\
--snp SNP \\
--a1 ALLELE1 \\
--a2 ALLELE0 \\
--signed-sumstats BETA,0 \\
--p P_BOLT_LMM_INF \\
--out ',loc_out_stem,'

echo "Running ldsc h2 estimation"
',dir_ldsc,'ldsc.py \\
--h2 ',paste0(loc_out_stem,".sumstats.gz"),' \\
--ref-ld-chr ',dir_LDscore,'LDscore. \\
--w-ld-chr ',dir_LDscore,'LDscore. \\
--out ',paste0(loc_out_stem,".h2"),'
                
', sep='', file=loc_script)
}


# then run inside script directory:
# for file in *.sh; do sbatch "$file";done
