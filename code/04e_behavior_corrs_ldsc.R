# computes phenotypic correlation between individual behaviors
# as well genetic correlations using LDsc

# Libraries and directories ####
library(tidyverse)
library(data.table)
library(polycor)
source('code/00_paths.R')

dir_ldsc_bvrs <- paste0(dir_scratch,'ldsc/behaviors/')
dir_ldsc <- '~/group_nuno/ldsc/' # ldsc software
dir_LDscore <- '~/group_nuno/LDscore/' # LD scores
dir.create(dir_ldsc_bvrs, showWarnings = FALSE)
dir_ldsc_pairs <- paste0(dir_ldsc_bvrs,'pairs/')
dir.create(dir_ldsc_pairs, showWarnings = FALSE)

fields <- as_tibble(fread(paste0(dir_scratch,'fields_tbl.txt')))
BRS_coeffs <- as_tibble(fread(paste0(dir_scratch,'BRS_models/BRS_coefficients.txt'))) %>%
  filter(term %in% fields$term[fields$use_type == 'behavior'],
         `beta-BRS-ALL-cov_bvr` != 0)
bvrs <- BRS_coeffs$term
pheno <- as_tibble(fread(paste0(dir_scratch,'cohorts/BRS_cohort_ALL.txt')))

# phenotypic correlations ####
col_BRS <- 'BRS-ALL-cov_bvr'
pheno_bvrs <- pheno %>% select(all_of(c(col_BRS,BRS_coeffs$term))) %>% as.data.frame()
# change columns to factors if categorical, booleans if binary
for (i in 1:length(colnames(pheno_bvrs))) {
  term <- colnames(pheno_bvrs)[i]
  print(paste(i,term))
  if (term == col_BRS) {next}
  data_type <- fields$data_type[fields$term == term]
  
  if (data_type == 'catord') { # ordered categorical variable
    # contains imputed values, which for these purposes are rounded to nearest whole integer
    pheno_bvrs[[term]] <- round(pheno_bvrs[[term]])
    # converts to factor, with levels in ascending order
    pheno_bvrs[[term]] <- factor(pheno_bvrs[[term]],
                                 levels = pheno_bvrs[[term]] %>% table() %>% names())
  } else if (data_type %in% c('catunord','binary')) {
    # contains imputed values, which for these purposes are rounded to nearest whole integer
    # converted to boolean
    pheno_bvrs[[term]] <- round(pheno_bvrs[[term]]) %>% as.logical()
  }
}
  

# gets correlation matrix for mix of continuous and categorical variables
hetcor_matrix <- hetcor(pheno_bvrs) # takes a while to run!

# combines correlation, standard error, and p-value into one table
cor_tibble <- as_tibble(as.data.frame(hetcor_matrix$correlations), rownames = "term1") %>%
  pivot_longer( cols = -term1, names_to = "term2", values_to = "r"
  ) %>% left_join(
    as_tibble(as.data.frame(hetcor_matrix$std.errors), rownames = "term1") %>%
      pivot_longer( cols = -term1, names_to = "term2", values_to = "r_se"
      ), by = c('term1','term2')
  ) %>% left_join(
    as_tibble(as.data.frame(hetcor_matrix$tests), rownames = "term1") %>%
      pivot_longer( cols = -term1, names_to = "term2", values_to = "p"
      ), by = c('term1','term2')
  ) %>% filter(term1 != term2)
loc_out <- paste0(dir_scratch, 'general_results/behavior_corr_table.tsv')
fwrite(cor_tibble, file=loc_out, sep='\t')



# ldsc prep LMM files ####

## behaviors + BRS ####
ldsc_files <- tibble(
  term = bvrs,
  term_ = str_replace_all(term, '\\.','_'),
  filepath = paste0(dir_scratch,'LMM/LMM.',term_,'.bgen.txt')
) %>% add_row( # manually insert BRS in there
  term = 'BRS', term_='BRS',
  filepath = paste0(dir_scratch,'LMM/LMM.ALL.BRS-ALL-cov_bvr.bgen.txt')
) %>% add_row( # manually insert CRS in there
  term = 'CRS', term_='CRS',
  filepath = paste0(dir_scratch,'LMM/LMM.CRS-ALL.bgen.txt')
) %>%
  mutate(
    sbatch_oe = paste0('ldsc.',term_),
    file_sh = paste0(sbatch_oe, '.sh'),
    loc_script = paste0(dir_ldsc_bvrs, file_sh),
    loc_out_stem = paste0(dir_ldsc_bvrs,'ldsc.',term_),
    loc_LMM.out = str_replace( filepath, '.bgen.txt','.out'),
    N = as.numeric(NA)
  ) 

###
# use this to filter the table to the actual scripts you need to currently run
#ldsc_files <- ldsc_files %>% filter(term == 'BRS')
###

for (i in 1:nrow(ldsc_files)) {
  slice <- ldsc_files[i,]
  print(paste(i, slice$term))
  
  # extracts sample size from .out file
  N <- read_lines(slice$loc_LMM.out) %>%
    str_subset("Number of individuals used in analysis: Nused =") %>%
    tail(1) %>% str_extract("(?<=Nused = )\\d+") %>% as.numeric()
  ldsc_files$N[i] <- N
  
  script <- cat('#!/bin/sh
#SBATCH -c 4
#SBATCH -t 0-00:10
#SBATCH -p short
#SBATCH --mem=4G
#SBATCH -o ',slice$sbatch_oe,'.out
#SBATCH -e ',slice$sbatch_oe,'.err

cd ',dir_ldsc_bvrs,'

module load miniconda3
conda env create --file ',dir_ldsc,'environment.yml
source activate ldsc

echo "Formatting summary file for ldsc"
',dir_ldsc,'munge_sumstats.py \\
--sumstats ',slice$filepath,' \\
--N ',N,' \\
--snp SNP \\
--a1 ALLELE1 \\
--a2 ALLELE0 \\
--signed-sumstats BETA,0 \\
--p P_BOLT_LMM_INF \\
--out ',slice$loc_out_stem,'

', sep='', file=slice$loc_script)
}



# then run the following bash command inside script directory:
cat('for file in',paste0(ldsc_files$file_sh, collapse=' '),'; do sbatch "$file";done')


## T2D ####

### prep summstats ####
# ok I'm breaking my whole rule about only having to run lower-numbered
# scripts, but I need you to run the bash code:
# 04f_downloading_T2DGWAS.sh
dir_T2DGWAS <- paste0(dir_scratch,'T2D/')
T2DGWAS_raw <- as_tibble(fread(paste0(dir_T2DGWAS,'All_Metal_LDSC-CORR_Neff.v2.txt')))
hg19_bed <- as_tibble(fread(paste0(dir_T2DGWAS,'rsID_coordinates_hg19.bed'))) %>%
  select(Chromsome=V1, Position=V3,SNP=V4) %>%
  mutate(Chromsome = str_replace(Chromsome,'chr','') %>% as.numeric())
# need an example of hg38 LMM
hg38_LMM <- as_tibble(fread(paste0(dir_scratch,'LMM/LMM.ALL.BRS-ALL-cov_bvr.bgen.txt'))) %>%
  select(Position_hg38=BP, SNP)

T2DGWAS_hg38 <- T2DGWAS_raw %>%
  left_join(hg19_bed, by=c('Chromsome','Position')) %>%
  left_join(hg38_LMM, by='SNP')

T2DGWAS <- T2DGWAS_hg38 %>%
  filter(!is.na(SNP)) %>%
  select(SNP, A1 = EffectAllele, A2 = NonEffectAllele,
         SIGNED_SUMSTATS = Beta, P=Pval, N_CAS_COL = Ncases, N_CON_COL = Ncontrols) %>%
  mutate(P = as.numeric(P))

# write hg38 T2DGWAS file to system (ready for ldsc's munge_sumstats)
loc_T2DGWAS <- paste0(dir_T2DGWAS,'T2DGWAS.txt')
fwrite(T2DGWAS, loc_T2DGWAS, sep='\t') 

### ldsc munge ####

# writes ldsc munge_sumstats script
script <- cat('#!/bin/sh
#SBATCH -c 4
#SBATCH -t 0-00:10
#SBATCH -p short
#SBATCH --mem=4G
#SBATCH -o ldsc.T2D.out
#SBATCH -e ldsc.T2D.err

cd ',dir_ldsc_bvrs,'

module load miniconda3
conda env create --file ',dir_ldsc,'environment.yml
source activate ldsc

echo "Formatting summary file for ldsc"
',dir_ldsc,'munge_sumstats.py \\
--sumstats ',loc_T2DGWAS,' \\
--N-cas-col N_CAS_COL \\
--N-con-col N_CON_COL \\
--snp SNP \\
--a1 A1 \\
--a2 A2 \\
--signed-sumstats SIGNED_SUMSTATS,0 \\
--p P \\
--merge-alleles ',paste0(dir_ldsc_bvrs,'ldsc.BRS.sumstats.gz'),' \\
--chunksize 1000000 \\
--out ',paste0(dir_ldsc_bvrs,'ldsc.T2D'),'

', sep='', file=paste0(dir_ldsc_bvrs,'ldsc.T2D.sh'))

# and submit that job from inside dir_ldsc_bvrs: sbatch ldsc.T2D.sh



## MAGIC HOMA-B and HOMA-IR ####

### prep summary stats ####

### ldsc munge ####
# from README pdf in MAGIC website
N_cases <- 51750
N_control <- 58074

for (HOMA in HOMAs) {
  # writes ldsc munge_sumstats script
  loc_HOMA <- paste0(dir_T2DGWAS, 'MAGIC_',HOMA,'_raw.txt')
  script <- cat('#!/bin/sh
#SBATCH -c 4
#SBATCH -t 0-00:10
#SBATCH -p short
#SBATCH --mem=4G
#SBATCH -o ldsc.',HOMA,'.out
#SBATCH -e ldsc.',HOMA,'.err

cd ',dir_ldsc_bvrs,'

module load miniconda3
conda env create --file ',dir_ldsc,'environment.yml
source activate ldsc

echo "Formatting summary file for ldsc"
',dir_ldsc,'munge_sumstats.py \\
--sumstats ',loc_HOMA,' \\
--N-cas ',N_cases,' \\
--N-con ',N_control,' \\
--snp Snp \\
--a1 effect_allele \\
--a2 other_allele \\
--signed-sumstats MainEffects,0 \\
--p MainP \\
--merge-alleles ',paste0(dir_ldsc_bvrs,'ldsc.BRS.sumstats.gz'),' \\
--chunksize 1000000 \\
--out ',paste0(dir_ldsc_bvrs,'ldsc.',HOMA),'

', sep='', file=paste0(dir_ldsc_bvrs,'ldsc.',HOMA,'.sh'))
}
# and submit both jobs from inside dir_ldsc_bvrs: sbatch ldsc.HOMAB.sh ldsc.HOMAIR.sh



# ldsc gencorrs ####

## behaviors x behaviors ####
# (plus HOMAB, HOMAIR, CRS, BRS, T2D, )

# makes table of all gencorr jobs
# each gencorr computation takes very little time (~60 seconds)
# but there are a lot of them: K*(K-1)/2 computations (for K behaviors)
# so I split the jobs across K-1 jobs, so each job is performing K/2 computations
bvrs2 <- c(HOMAs,'CRS','T2D','BRS',bvrs)
K <- length(bvrs2)
ldsc_rgs <- combn(str_replace_all(bvrs2,'\\.','_'), 2) %>%
  t() %>% as_tibble() %>% rename(term1_=V1, term2_=V2) %>% mutate(
    filepath1 = paste0(dir_ldsc_bvrs,'ldsc.',term1_,'.sumstats.gz'),
    filepath2 = paste0(dir_ldsc_bvrs,'ldsc.',term2_,'.sumstats.gz'),
    loc_out_stem = paste0('ldsc.rg.',term1_,'.',term2_),
    k = ceiling( row_number() / (K/2) )
  )

max_mins_per_rg = 2
ldsc_jobs <- tibble(
  k = 1:(K-1),
  job_i_start = sapply(1:(K-1), function(k) min(which(ldsc_rgs$k == k))),
  job_i_stop = sapply(1:(K-1), function(k) max(which(ldsc_rgs$k == k))),
  sbatch_oe = paste0('ldsc.rg.',k),
  file_sh = paste0(sbatch_oe ,'.sh'),
  loc_script = paste0(dir_ldsc_pairs, file_sh),
  computations = table( ldsc_rgs$k ) %>% c(),
  compute_mins = round(computations * max_mins_per_rg) # keep under 60
  )


for (k in 1:(K-1)) {
  print(k)
  slice_k <- ldsc_jobs[k,]
  
  script <- cat('#!/bin/sh
#SBATCH -c 4
#SBATCH -t 0-00:',slice_k$compute_mins,'
#SBATCH -p short
#SBATCH --mem=4G
#SBATCH -o ',slice_k$sbatch_oe,'.out
#SBATCH -e ',slice_k$sbatch_oe,'.err

cd ',dir_ldsc_pairs,'

module load miniconda3
#conda env create --file ',dir_ldsc,'environment.yml
source activate ldsc
', sep='', file=slice_k$loc_script)
  
  for (i in slice_k$job_i_start:slice_k$job_i_stop) {
    slice_i <- ldsc_rgs[i,]
    script_rg <- cat('
echo "Running LDscore genetic correlation"
echo ',slice_i$term1_,',',slice_i$term2_,'
',dir_ldsc,'ldsc.py \\
--rg ',slice_i$filepath1,',',slice_i$filepath2,' \\
--ref-ld-chr ',dir_LDscore,'LDscore. \\
--w-ld-chr ',dir_LDscore,'LDscore. \\
--out ',slice_i$loc_out_stem,'
', sep='', file=slice_k$loc_script, append=TRUE)
  }
}

# then run the following bash command inside script directory:
cat('for file in',paste0(ldsc_jobs$file_sh, collapse=' '),'; do sbatch "$file";done')