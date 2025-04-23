# makes table summarizing analyses for behaviors + BRS

# loads directories and packages ####
library(tidyverse)
library(data.table)
source('code/00_paths.R')

# relevant files/variables
ROC_tbl <- as_tibble(fread(paste0(dir_scratch,'general_results/ROC_tbl.tsv')))
shortnames <- as_tibble(fread(paste0(dir_repo,'input_data/term_shortnames.tsv')))
gencorr_REML <- as_tibble(fread(paste0(dir_scratch,'general_results/gencorr_REML.tsv')))
ldsc_rg <- as_tibble(fread(paste0(dir_scratch,'general_results/behavior_rg_table.tsv')))
GWAS_summary_tbl <- as_tibble(fread(paste0(dir_scratch,'general_results/GWAS_summary_tbl.txt')))
gene_assoc_n <- as_tibble(fread(paste0(dir_scratch,'general_results/gene_associations_n.tsv')))

fields <- as_tibble(fread(paste0(dir_scratch,'fields_tbl.txt')))
col_bvrs <- fields$term[fields$use_type == 'behavior']
col_BRS <- 'BRS-ALL-cov_bvr' # main BRS term
ROC_tbl$term[ROC_tbl$term == 'BRS'] <- col_BRS

# loads original cvglm BRS model and coefficients
load(paste0(dir_scratch,'BRS_models/cv_glm_models.RData'))
BRS_cvglm <- models[[str_remove(col_BRS, 'BRS-')]]
BRS_coeffs <- as_tibble(fread(paste0(dir_scratch,'BRS_models/BRS_coefficients.txt'))) %>%
  select(term, traitname, beta = all_of(paste0('beta-',col_BRS))) %>%
  drop_na() %>% 
  mutate(factor_SD = BRS_cvglm$factor_SDs,
         beta_norm = beta * factor_SD,
         beta_norm_abs = abs(beta_norm),
         effect_sign = sign(beta_norm),
         hazard_ratio = exp(beta_norm)) %>%
  filter(beta != 0)

# constructs summary table (behaviors + BRS)
sumtbl <- tibble(term = c(col_BRS, col_bvrs[col_bvrs %in% BRS_coeffs$term])) %>%
  # adds shortened names
  left_join(shortnames %>%
              mutate(category = ifelse(category=='',NA,category)), by='term') %>%
  # adds data and use type
  left_join(fields %>% select(term, data_type, use_type), by='term') %>%
  # adds BRS coefficient (set to 1 for BRS itself)
  left_join(BRS_coeffs %>% select(term, beta, beta_norm) %>%
              add_row(term=col_BRS, beta=1, beta_norm=1), by='term') %>%
  # adds AUC for T2D
  left_join(ROC_tbl, by='term') %>%
  # adds h2 estimate and rg w/ CRS-ALL estimate (BOLT-REML)
  left_join(gencorr_REML %>% 
              mutate(term = str_replace_all(trait2,'_','\\.')) %>%
              filter(term %in% BRS_coeffs$term, trait1 == 'CRS-ALL') %>%
              mutate(h2=trait2_h2g, h2_se = trait2_h2g_err) %>%
              add_row(
                gencorr_REML %>% filter(trait1==col_BRS, trait2=='CRS-ALL') %>%
                  mutate(term=trait1, h2=trait1_h2g, h2_se = trait1_h2g_err)) %>%
              select(term, h2, h2_se, rg_CRS = rg, rg_se_CRS = rg_err), by='term') %>%
  # adds ldsc: h2 estimate, intercept
  left_join(ldsc_rg %>% filter(term1 == 'T2D') %>%
              select(term=term2, h2_ldsc = h2g2, h2_se_ldsc = h2g_se2,
                     ldsc_intercept=intercept2) %>%
              mutate(term = ifelse(term == 'BRS', col_BRS,
                                   str_replace(term,'_','.'))),
            by='term') %>%
  # adds ldsc rg w/ CRS, HOMA-B, HOMA-IR
  left_join(ldsc_rg %>%
              filter(term1 %in% c('T2D','CRS', 'HOMAB', 'HOMAIR')) %>%
              mutate(term1 = ifelse(term1=='CRS','CRS_ldsc',term1)) %>%
              select(term=term2, term1, rg, rg_se) %>%
              pivot_wider(names_from = term1,
                          names_glue = '{.value}_{term1}',
                          values_from = c(rg, rg_se)) %>%
              mutate(term = ifelse(term == 'BRS', col_BRS,
                                   str_replace(term,'_','.'))),
            by='term') %>%
  # adds GWAS summary data
  left_join(GWAS_summary_tbl %>%
              mutate(term = ifelse(trait == 'ALL.BRS-ALL-cov_bvr', col_BRS, trait)) %>%
              select(-trait), by='term') %>%
  # adds # of genes for GWAS signals
  left_join(gene_assoc_n %>% rename(n_sig_genes = n), by='shortname') %>% 
  mutate(term = ifelse(term == col_BRS, 'BRS', term))



# writes to folder
fwrite(sumtbl, paste0(dir_scratch,'general_results/sumtbl.tsv'), sep='\t')
