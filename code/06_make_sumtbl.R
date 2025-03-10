# makes table summarizing analyses for behaviors + BRS

# loads directories and packages ####
library(tidyverse)
library(data.table)
source('code/00_paths.R')

# relevant files/variables
gencorr_REML <- as_tibble(fread('scratch/general_results/gencorr_REML.tsv'))
ROC_tbl <- as_tibble(fread('scratch/general_results/ROC_tbl.tsv'))
shortnames <- as_tibble(fread('input_data/term_shortnames.tsv'))
fields <- as_tibble(fread('scratch/fields_tbl.txt'))
col_bvrs <- fields$term[fields$use_type == 'behavior']
col_BRS <- 'BRS-ALL-cov_bvr' # main BRS term
ROC_tbl$term[ROC_tbl$term == 'BRS'] <- col_BRS

# loads original cvglm BRS model and coefficients
load('scratch/BRS_models/cv_glm_models.RData')
BRS_cvglm <- models[[str_remove(col_BRS, 'BRS-')]]
BRS_coeffs <- as_tibble(fread('scratch/BRS_models/BRS_coefficients.txt')) %>%
  select(term, traitname, beta = all_of(paste0('beta-',col_BRS))) %>%
  drop_na() %>% 
  mutate(factor_SD = BRS_cvglm$factor_SDs,
         beta_norm = beta * factor_SD,
         beta_norm_abs = abs(beta_norm),
         effect_sign = sign(beta_norm)) %>%
  filter(beta != 0)

# constructs summary table (behaviors + BRS)
sumtbl <- tibble(term = c(col_BRS, col_bvrs[col_bvrs %in% BRS_coeffs$term])) %>%
  # adds shortened names
  left_join(shortnames, by='term') %>%
  # adds data and use type
  left_join(fields %>% select(term, data_type, use_type), by='term') %>%
  # adds BRS coefficient (set to 1 for BRS itself)
  left_join(BRS_coeffs %>% select(term, beta, beta_norm) %>%
              add_row(term=col_BRS, beta=1, beta_norm=1), by='term') %>%
  # adds AUC for T2D
  left_join(ROC_tbl, by='term') %>%
  # adds h2 estimate and rg w/ CRS-ALL estimate
  left_join(gencorr_REML %>% 
              mutate(term = str_replace_all(trait2,'_','\\.')) %>%
              filter(term %in% BRS_coeffs$term, trait1 == 'CRS-ALL') %>%
              mutate(h2=trait2_h2g, h2_err = trait2_h2g_err) %>%
              add_row(
                gencorr_REML %>% filter(trait1==col_BRS, trait2=='CRS-ALL') %>%
                  mutate(term=trait1, h2=trait1_h2g, h2_err = trait1_h2g_err)) %>%
              select(term, h2, h2_err, rg_CRS = rg, rg_err_CRS = rg_err), by='term')

# writes to folder
fwrite(sumtbl, 'scratch/general_results/sumtbl.tsv', sep='\t')
