# makes final tables for manuscript

# loads directories and packages ####
library(tidyverse)
library(data.table)
source('code/00_paths.R')
source('code/00_plotting.R')
col_BRS <- 'BRS-ALL-cov_bvr' # main BRS term

sumtbl <- as_tibble(fread(paste0(dir_scratch,'general_results/sumtbl.tsv')))
fields <- as_tibble(fread(paste0(dir_scratch,'fields_tbl.txt')))

# all_behaviors.csv ####
all_behaviors <- as_tibble(fread(paste0(dir_repo,'input_data/UKBB_exposures.txt'))) %>%
  filter(Exposure_Class == 'agency') %>%
  select(-Exposure_Class)
loc_out <- paste0(dir_tbls,'all_behaviors.csv')
fwrite(all_behaviors, loc_out, sep=',')

# BRS_coefficients.csv ####
load(paste0(dir_scratch,'BRS_models/cv_glm_models.RData'))
BRS_cvglm <- models[[str_remove(col_BRS, 'BRS-')]]
BRS_coefficients <- as_tibble(fread(paste0(dir_scratch,'BRS_models/BRS_coefficients.txt'))) %>%
  left_join(sumtbl %>% select(term, shortname), by='term') %>%
  select(term, shortname, beta = all_of(paste0('beta-',col_BRS))) %>%
  filter(!is.na(beta)) %>%
  mutate(factor_SD = BRS_cvglm$factor_SDs,
         beta_per_SD = beta * factor_SD,) %>%
  filter(beta != 0) %>% select(-factor_SD)
loc_out <- paste0(dir_tbls,'BRS_coefficients.csv')
fwrite(BRS_coefficients, loc_out, sep=',')


# behavior_summary_results.csv ####
behavior_summary_results <- sumtbl %>% 
  left_join(fields %>% select(term, fieldID, fieldValue = value), by='term') %>%
  select(term, fieldID, fieldValue, everything(),
         -traitname, -beta, -beta_norm, -use_type)
loc_out <- paste0(dir_tbls,'behavior_summary_results.csv')
fwrite(behavior_summary_results, loc_out, sep=',')

# model_AUCs.csv ####
ROC_tbl <- as_tibble(fread(paste0(dir_scratch, 'general_results/ROC_tbl.tsv')))
shortnames <- as_tibble(fread(paste0(dir_repo,'input_data/term_shortnames.tsv'))) %>%
  distinct() %>%
  add_row(term='cov', shortname='Covariate-only')
shortnames$term[shortnames$term == 'PRS_T2D'] <- 'PRS'

model_AUCs <- ROC_tbl %>%
  left_join(shortnames %>% select(term, shortname), by='term') %>%
  mutate(shortname = ifelse(is.na(shortname), str_replace_all(term, '\\+',' + '), shortname)) %>%
  select(term, shortname, AUC, AUC_low, AUC_upp, AUC_se)
loc_out <- paste0(dir_tbls,'model_AUCs.csv')
fwrite(model_AUCs, loc_out, sep=',')

