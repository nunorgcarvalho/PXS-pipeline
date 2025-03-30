# reads ldsc h2 data

library(tidyverse)
library(data.table)
source('code/00_paths.R')

dir_ldsc_out <- paste0(dir_scratch,'ldsc/')
ldsc_files <- list.files(paste0(dir_ldsc_out), pattern='*.h2.log')

ldsc_h2 <- tibble(trait='',N_SNPs=0, h2=0,h2_se=0,lambda_GC=0,intercept=0,intercept_se=0)[0,]

for (ldsc_file in ldsc_files) {
  filepath <- paste0(dir_ldsc_out, ldsc_file)
  trait <- str_replace(str_sub(ldsc_file, 6),'.h2.log','')
  
  h2_log <- read_lines(filepath)
  
  # extracts statistics
  ldsc_h2 <- ldsc_h2 %>% add_row(
    trait=trait,
    N_SNPs=(h2_log %>% str_subset("with regression SNP LD,") %>% str_split(' '))[[1]][7] %>% as.numeric(),
    h2=(h2_log %>% str_subset("scale h2: ") %>% str_split(' '))[[1]][5] %>% as.numeric(),
    h2_se=(h2_log %>% str_subset("scale h2: ") %>% str_split(' '))[[1]][6] %>% str_sub(2,-2) %>% as.numeric(),
    lambda_GC=(h2_log %>% str_subset("Lambda GC: ") %>% str_split(' '))[[1]][3] %>% as.numeric(),
    intercept=(h2_log %>% str_subset("Intercept: ") %>% str_split(' '))[[1]][2] %>% as.numeric(),
    intercept_se=(h2_log %>% str_subset("Intercept: ") %>% str_split(' '))[[1]][3] %>% str_sub(2,-2) %>% as.numeric()
  )
}


# Visualizing table
ldsc_h2_clean <- ldsc_h2 %>%
  mutate(h2_cohort = sapply(str_split(trait, "\\."), `[`, 1),
         BRS = sapply(str_split(trait, "\\."), `[`, 2),
         training_cohort = sapply(str_split(trait, "-"), `[`, 2),
         model_factors = sapply(str_split(trait, "-"), `[`, 3),
         model_label = str_replace_all(model_factors, '_',' + ')) %>%
  select(h2_cohort, training_cohort, model_label, h2, h2_se, lambda_GC, intercept, intercept_se)





# rg pairs of behaviors ####
dir_pairs <- paste0(dir_ldsc_out,'behaviors/pairs/')
ldsc_files <- list.files(dir_pairs, pattern='*.log')

ldsc_rg <- tibble(
  file = ldsc_files,
  term1 = map_chr( str_split( ldsc_files, '\\.'), 3),
  term2 = map_chr( str_split( ldsc_files, '\\.'), 4),
  h2g1=NA, h2g_se1=NA, h2g2=NA, h2g_se2=NA, rg=NA, rg_se=NA,
  intercept1=NA, intercept_se1=NA, intercept2=NA, intercept_se2=NA)

for (i in 1:nrow(ldsc_rg)) {
  if (i %% 100 == 0) {print(i)}
  filepath <- paste0(dir_pairs, ldsc_rg$file[i])
  rg_log <- read_lines(filepath)
  
  # extracts info from ldsc log file
  h2gs <- map_chr( rg_log %>% str_subset('Observed scale h2') %>% str_split(' '), 5 ) %>% as.numeric()
  h2g_ses <- map_chr( rg_log %>% str_subset('Observed scale h2') %>% str_split(' '), 6 ) %>% str_sub(2,-2) %>% as.numeric()
  intercepts <- map_chr( rg_log %>% str_subset('Intercept: ') %>% str_split(' '), 2 ) %>% as.numeric()
  intercept_ses <- map_chr( rg_log %>% str_subset('Intercept: ') %>% str_split(' '), 2 ) %>% as.numeric()
  rg <- map_chr( rg_log %>% str_subset('Genetic Correlation: ') %>% str_split(' '), 3 ) %>% as.numeric()
  rg_se <- map_chr( rg_log %>% str_subset('Genetic Correlation: ') %>% str_split(' '), 4 ) %>% str_sub(2,-2) %>% as.numeric()
  # adds to table
  ldsc_rg[i,4:13] <- list(h2gs[1],h2g_ses[1],h2gs[2],h2g_ses[2],rg,rg_se,
                          intercepts[1],intercept_ses[1], intercepts[2],intercept_ses[2])
  
}
# saves to system
loc_out <- paste0(dir_scratch, 'general_results/behavior_rg_table.tsv')
fwrite(ldsc_rg, file=loc_out, sep='\t')
