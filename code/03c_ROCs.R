# script for comparing (AU)ROC for behaviors vs BRS

# loads directories and packages ####
library(tidyverse)
library(data.table)
library(pROC)
source('code/00_paths.R')

# general files/variables ####
col_BRS <- 'BRS-ALL-cov_bvr' # main BRS term
shortnames <- as_tibble(fread('input_data/term_shortnames.tsv'))
fields <- as_tibble(fread('scratch/fields_tbl.txt'))
col_covs_ALL <- fields$term[fields$use_type == 'covar']
col_bvrs_ALL <- fields$term[fields$use_type == 'behavior']
col_CRFs <- fields$term[fields$use_type == 'CRF']

# reads phenotype and relevant columns
pheno <- as_tibble(fread('scratch/cohorts/BRS_cohort_ALL.txt'))
BRS_coeffs <- as_tibble(fread('scratch/BRS_models/BRS_coefficients.txt')) %>%
  select(term, traitname, beta = all_of(paste0('beta-',col_BRS))) %>%
  drop_na() %>% filter(beta != 0)
col_covs <- BRS_coeffs$term[BRS_coeffs$term %in% col_covs_ALL]
col_bvrs <- BRS_coeffs$term[BRS_coeffs$term %in% col_bvrs_ALL]

# i end up only using the sex, age, PC1,11-12, because these were the only
# covariates that came up as significant. saves some compute time


# defines training and testing subsets
training_tbl <- pheno %>% filter(sample_group == 'A')
testing_tbl <- pheno %>% filter(sample_group == 'B')
test_y <- testing_tbl$T2D_onset # for ROC later

# reused from 02a_compute_BRS.R
get_cvglm_obj <- function(training_tbl, cols_use) {
  training_tbl <- training_tbl[,c('T2D_onset_days','T2D_onset',cols_use)] %>% drop_na()
  glm_y <- survival::Surv(training_tbl$T2D_onset_days,training_tbl$T2D_onset)
  glm_x <- as.matrix( training_tbl[,cols_use] )
  
  cv_glm1 <- glmnet::cv.glmnet(glm_x, glm_y, family='cox', alpha=1, nfolds=10,
                               type.measure = 'C')
  return(cv_glm1)
}

model_factors <- list(
  'cov' = c(col_covs),
  'BRS' = c(col_covs, col_bvrs_ALL),
  'CRS' = c(col_covs, col_CRFs[-4]), # excludes HbA1c
  'PRS' = c(col_covs, 'PRS_T2D'),
  'BRS+CRS' = c(col_covs, col_bvrs_ALL, col_CRFs[-4]),
  'BRS+PRS' = c(col_covs, col_bvrs_ALL, 'PRS_T2D'),
  'CRS+PRS' = c(col_covs, col_CRFs[-4], 'PRS_T2D'),
  'BRS+CRS+PRS' = c(col_covs, col_bvrs_ALL, col_CRFs[-4], 'PRS_T2D')
)

ROC_tbl <- tibble(term = c('cov',col_bvrs,'BRS',col_CRFs,'CRS','PRS',
                           'BRS+CRS','BRS+PRS','CRS+PRS','BRS+CRS+PRS'))

cvglm_list <- list()
ROC_list <- list()

set.seed(2025)
for (i in 1:nrow(ROC_tbl)) {
  term <- ROC_tbl$term[i]
  print(paste(i, term))
  
  if (term %in% names(model_factors)) {
    cols_use <- model_factors[[term]]
  } else {
    cols_use <- c(col_covs, term)
  }
  
  if (i >= 1) {
    cvglm_i <- get_cvglm_obj(training_tbl, cols_use)
    cvglm_list[[term]] <- cvglm_i
  } else {
    cvglm_i <- cvglm_list[[term]]
  }
  
  betas <- coef(cvglm_i, s='lambda.1se') %>% as.numeric()
  factors_matrix <- as.matrix( testing_tbl[,cols_use] )
  
  pred_raw <- (factors_matrix %*% betas)[,1]
  pred_vec <- scale(pred_raw)[,1]
  
  roc1 <- roc(test_y, pred_vec, ci=TRUE)
  ROC_list[[term]] <- roc1
  
}

#save(cvglm_list, file=paste0(dir_scratch,'BRS_models/cvglm_list_terms.RData'))
#save(ROC_list, file=paste0(dir_scratch,'BRS_models/ROC_list_terms.RData'))

#load(paste0(dir_scratch,'BRS_models/cvglm_list_terms.RData'))
#load(paste0(dir_scratch,'BRS_models/ROC_list_terms.RData'))



# adds AUC information to general table
ROC_tbl$AUC <- map_dbl(ROC_list, 'auc') # extracts AUC
ROC_tbl$AUC_low <- map_dbl(ROC_list, ~ as.numeric(.$ci)[1] )
ROC_tbl$AUC_upp <- map_dbl(ROC_list, ~ as.numeric(.$ci)[3] )
CI95_z <- qnorm(1 - (1 - 0.95)/2)
ROC_tbl$AUC_se <- (ROC_tbl$AUC_upp - ROC_tbl$AUC_low) / (2*CI95_z)

# saves to system
fwrite(ROC_tbl, file=paste0(dir_scratch, 'general_results/ROC_tbl.tsv'))

# View top-performing models
ROC_tbl %>% left_join(shortnames, by='term') %>%
  mutate(shortname = ifelse(is.na(shortname), term, shortname)) %>%
  select(term, shortname, AUC, AUC_low, AUC_upp) %>%
  arrange(-AUC) %>% View()


# median ROC
median( (ROC_tbl %>% filter(term %in% col_bvrs))$AUC)



# makes long table of ROC points (sensitivity, specificity, threshold) for each term
ROC_long <- tibble(term='',sensitivity=0,specificity=0,threshold=0)[0,]
max_points <- 200 # how many points are plotted per ROC
set.seed(2025)
for (i in 1:length(ROC_list)) {
  term <- names(ROC_list)[i]
  ROC_obj <- ROC_list[[term]]
  keep_i <- 1:length(ROC_obj$thresholds)
  # randomly samples `max_threshold` # of ROC points, including first & last
  if (length(ROC_obj$thresholds) > max_points) {
    keep_i <- c(1, length(ROC_obj$thresholds),
                sample(1:length(ROC_obj$thresholds), max_points - 2)) %>%
      sort()
  }
  ROC_long <- ROC_long %>% add_row(
    term=term, sensitivity=ROC_obj$sensitivities[keep_i],
    specificity=ROC_obj$specificities[keep_i],
    threshold = ROC_obj$thresholds[keep_i]
  )
}
# calculates Youden J statistic, and max per term
ROC_long$youden_J <- ROC_long$sensitivity + ROC_long$specificity - 1
ROC_youden <- ROC_long %>% group_by(term) %>%
  summarize(max_youden_J = max(youden_J))


# BRS vs bvrs ####
col_highlight <- c('BRS', 'cov')
# makes table for BRS vs behavior plotting purposes
ROC_long_bvr <- ROC_long %>% filter(term %in% c(col_bvrs, col_highlight) ) %>%
  mutate(color = ifelse(term == 'BRS','red',
                        ifelse(term == 'cov','black','gray60')),
         size = ifelse(term %in% col_highlight, 0.5, 0.1)) %>%
  arrange(term %in% col_highlight) # easier plotting later
ROC_youden_bvr <- ROC_long_bvr %>%
  left_join(ROC_youden, by='term') %>%
  filter(youden_J == max_youden_J) %>%
  mutate(size = ifelse(term=='BRS', 2, 1))

# plot
ggplot(ROC_long_bvr, aes(x = 1 - specificity, y = sensitivity,
                         color=color, size=size)) +
  geom_line(data=ROC_long_bvr %>% filter(!(term %in% col_highlight)), aes(group=term)) +
  geom_line(data=ROC_long_bvr %>% filter(term %in% col_highlight), aes(group=term)) +
  geom_abline(slope=1, linetype='dashed') +
  geom_point(data=ROC_youden_bvr) +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  scale_color_identity() + scale_size_identity() +
  theme_bw() + theme(legend.position = 'none') +
  labs(x = '1 - Specificity',
       y = 'Sensitivity',
       title = 'The Behavioral Risk Score predicts Type 2 Diabetes better than individual behaviors',
       subtitle = paste0(
         'Youden point for each ROC is marked',
         '\nIndividual behaviors in gray. Behavioral Risk Score (BRS) in red. Covariate-only model in black.'))

# View 
shortnames$term[shortnames$term == col_BRS] <- 'BRS'
ROC_tbl %>% filter(term %in% c(col_bvrs, 'BRS') ) %>% arrange(-AUC) %>%
  left_join(shortnames, by='term') %>%
  select(term, shortname, starts_with('AUC'))


# CRS vs CRFs ####
col_highlight <- c('CRS', 'cov')
# makes table for CRS vs CRFs plotting purposes
ROC_long_CRF <- ROC_long %>% filter(term %in% c(col_CRFs[-4], col_highlight) ) %>%
  mutate(color = ifelse(term == 'CRS','red',
                        ifelse(term == 'cov','black','gray60')),
         size = ifelse(term %in% col_highlight, 0.5, 0.1)) %>%
  arrange(term %in% col_highlight) # easier plotting later
ROC_youden_CRF <- ROC_long_CRF %>%
  left_join(ROC_youden, by='term') %>%
  filter(youden_J == max_youden_J) %>%
  mutate(size = ifelse(term=='CRS', 2, 1))

# plot
ggplot(ROC_long_CRF, aes(x = 1 - specificity, y = sensitivity,
                         color=color, size=size)) +
  geom_line(data=ROC_long_CRF %>% filter(!(term %in% col_highlight)), aes(group=term)) +
  geom_line(data=ROC_long_CRF %>% filter(term %in% col_highlight), aes(group=term)) +
  geom_abline(slope=1, linetype='dashed') +
  geom_point(data=ROC_youden_CRF) +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  scale_color_identity() + scale_size_identity() +
  theme_bw() + theme(legend.position = 'none') +
  labs(x = '1 - Specificity',
       y = 'Sensitivity')



# risk score comparisons ####
ROC_tbl_scores <- ROC_tbl %>%
  filter(!(term %in% c(col_bvrs_ALL, col_covs_ALL, col_CRFs))) %>%
  mutate(shortname = ifelse(term == 'cov', 'covariate-only',
         str_replace_all(term, '\\+',' + '))) %>%
  mutate(shortname = factor(shortname, levels=shortname))

ggplot(ROC_tbl_scores, aes(x=fct_rev(shortname))) +
  geom_point(aes(y=AUC)) +
  geom_errorbar(aes(ymin = AUC_low, ymax = AUC_upp),
                width=0.5) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,1/20),
                     minor_breaks = seq(0,1,1/100)) +
  labs(x = 'Model',
       y = 'AUC (95% CI)') +
  coord_flip()

## DeLong test for comparing ROCs ####
roc.test(ROC_list[['BRS']], ROC_list[['PRS']])
roc.test(ROC_list[['BRS']], ROC_list[['CRS']])
roc.test(ROC_list[['BRS']], ROC_list[['f21001']])
roc.test(ROC_list[['BRS']], ROC_list[['BRS+PRS']])
roc.test(ROC_list[['CRS']], ROC_list[['BRS+CRS']])
roc.test(ROC_list[['CRS']], ROC_list[['CRS+PRS']])
roc.test(ROC_list[['BRS+CRS']], ROC_list[['CRS+PRS']])
roc.test(ROC_list[['BRS+CRS+PRS']], ROC_list[['CRS+PRS']])
roc.test(ROC_list[['BRS+CRS+PRS']], ROC_list[['f30750']])
