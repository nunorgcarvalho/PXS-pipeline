# script for comparing (AU)ROC for behaviors vs BRS

# loads directories and packages ####
library(tidyverse)
library(data.table)
library(pROC)
library(ggrepel)
library(ggpubr)
source('code/00_paths.R')
source('code/00_plotting.R')

# general files/variables ####
col_BRS <- 'BRS-ALL-cov_bvr' # main BRS term
shortnames <- as_tibble(fread(paste0(dir_repo,'input_data/term_shortnames.tsv')))
fields <- as_tibble(fread(paste0(dir_scratch,'fields_tbl.txt')))
col_covs_ALL <- fields$term[fields$use_type == 'covar']
col_bvrs_ALL <- fields$term[fields$use_type == 'behavior']
col_CRFs <- fields$term[fields$use_type == 'CRF']

# reads phenotype and relevant columns
BRS_coeffs <- as_tibble(fread(paste0(dir_scratch,'BRS_models/BRS_coefficients.txt'))) %>%
  select(term, traitname, beta = all_of(paste0('beta-',col_BRS))) %>%
  drop_na() %>% filter(beta != 0)
col_covs <- BRS_coeffs$term[BRS_coeffs$term %in% col_covs_ALL]
col_bvrs <- BRS_coeffs$term[BRS_coeffs$term %in% col_bvrs_ALL]
pheno <- as_tibble(fread(paste0(dir_scratch,'cohorts/BRS_cohort_ALL.txt')))

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
  cvglm_i <- get_cvglm_obj(training_tbl, cols_use)
  cvglm_list[[term]] <- cvglm_i
  
  # gets standard deviations of factors, for later rescaling if desired
  cvglm_list[[term]]$factor_SDs <- sapply(training_tbl[,cols_use], function(x) sd(x, na.rm=TRUE))
  
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
#fwrite(ROC_tbl, file=paste0(dir_scratch, 'general_results/ROC_tbl.tsv'))
#ROC_tbl <- as_tibble(fread(paste0(dir_scratch, 'general_results/ROC_tbl.tsv')))

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
  mutate(color_label = ifelse(term == 'BRS', 'Behavioral Risk Score',
                              ifelse(term == 'cov', 'Covariate-only', 'Behavior')),
         size = ifelse(term %in% col_highlight, 0.5, 0.1)) %>%
  arrange(term %in% col_highlight) # easier plotting later
ROC_youden_bvr <- ROC_long_bvr %>%
  left_join(ROC_youden, by='term') %>%
  filter(youden_J == max_youden_J) %>%
  mutate(size = ifelse(term=='BRS', 2, 1))

# ROC plot shared features
add_ROC_plot_elements <- function(gg, legend=TRUE) {
  gg <- gg +
    geom_abline(slope=1, linetype='dashed') +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    scale_size_identity() +
    labs(color = 'Model Type',
         x = '1 - Specificity',
         y = 'Sensitivity') +
    coord_fixed(expand = FALSE) +
    theme_bw() +
    theme(legend.position = c(1,0), legend.justification = c(1.05,-0.05),
          legend.box.background = element_rect(color = "black"),
          plot.margin = margin(r = 5, l=1, unit = "mm"))
  if (legend==FALSE) {
    gg <- gg + theme(legend.position = 'none')
  }
  return(gg)
}

# plot
gg <- ggplot(ROC_long_bvr, aes(x = 1 - specificity, y = sensitivity,
                         color=color_label, size=size)) +
  geom_line(data=ROC_long_bvr %>% filter(!(term %in% col_highlight)), aes(group=term)) +
  geom_line(data=ROC_long_bvr %>% filter(term %in% col_highlight), aes(group=term)) +
  geom_point(data=ROC_youden_bvr) +
  scale_color_manual(breaks = c('Behavioral Risk Score','Behavior','Covariate-only'),
                     values = c('red','gray60','black'))
add_ROC_plot_elements(gg)
  
loc_fig <- paste0(dir_figs,'ROC_bvrs_BRS')
ggsave(paste0(loc_fig,".png"), width=180, height=180, units="mm", dpi=300)

# View 
ROC_tbl %>% filter(term %in% c(col_bvrs, 'BRS') ) %>% arrange(-AUC) %>%
  left_join(shortnames, by='term') %>%
  select(term, shortname, starts_with('AUC'))


# CRS vs CRFs ####
col_highlight <- c('CRS', 'cov')
# makes table for CRS vs CRFs plotting purposes
ROC_long_CRF <- ROC_long %>% filter(term %in% c(col_CRFs[-4], col_highlight) ) %>%
  mutate(color_label = ifelse(term == 'CRS', 'Clinical Risk Score',
                              ifelse(term == 'cov', 'Covariate-only', 'Clinical Risk Factor')),
         size = ifelse(term %in% col_highlight, 0.5, 0.1)) %>%
  arrange(term %in% col_highlight) # easier plotting later
ROC_youden_CRF <- ROC_long_CRF %>%
  left_join(ROC_youden, by='term') %>%
  filter(youden_J == max_youden_J) %>%
  mutate(size = ifelse(term=='CRS', 2, 1)) %>%
  left_join(shortnames %>% select(term, shortname) %>%
              filter(term %in% col_CRFs), by='term')

# plot
gg <- ggplot(ROC_long_CRF, aes(x = 1 - specificity, y = sensitivity,
                         color=color_label, size=size)) +
  geom_line(data=ROC_long_CRF %>% filter(!(term %in% col_highlight)), aes(group=term)) +
  geom_line(data=ROC_long_CRF %>% filter(term %in% col_highlight), aes(group=term)) +
  geom_point(data=ROC_youden_CRF) +
  geom_text_repel(data=ROC_youden_CRF, aes(label=shortname),
                  color='gray60', size=2) +
  scale_color_manual(breaks = c('Clinical Risk Score','Clinical Risk Factor','Covariate-only'),
                     values = c('red','gray60','black'))
add_ROC_plot_elements(gg)

loc_fig <- paste0(dir_figs,'ROC_CRF_CRS')
ggsave(paste0(loc_fig,".png"), width=180, height=180, units="mm", dpi=300)




# risk score comparisons ####

## AUC comparison ####
ROC_tbl_scores <- ROC_tbl %>%
  filter(!(term %in% c(col_bvrs_ALL, col_covs_ALL, col_CRFs))) %>%
  mutate(shortname = ifelse(term == 'cov', 'Covariate-only',
                     ifelse(term == 'BRS', 'Behavioral Risk Score (BRS)',
                     ifelse(term == 'CRS', 'Clinical Risk Score (CRS)',
                     ifelse(term == 'PRS', 'Polygenic Risk Score (PRS)',
                            str_replace_all(term, '\\+',' + ')))))) %>%
  mutate(shortname = factor(shortname, levels=shortname))

gg_mm1 <- ggplot(ROC_tbl_scores, aes(x=fct_rev(shortname), color=shortname)) +
  geom_point(aes(y=AUC)) +
  geom_errorbar(aes(ymin = AUC_low, ymax = AUC_upp),
                width=0.3) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,1/20), minor_breaks = seq(0,1,1/100)) +
  labs(x = 'Model',
       y = 'AUC for Type 2 Diabetes (95% CI)') +
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width=15)) +
  theme(legend.position='none',
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))
loc_fig <- paste0(dir_figs,"main_model_AUCs")
ggsave(paste0(loc_fig,".png"), gg_mm1, width=120, height=180, units="mm", dpi=300)

## ROC comparison ####
# makes table for CRS vs CRFs plotting purposes
ROC_long_RS <- ROC_long %>% filter(term %in% ROC_tbl_scores$term ) %>%
  mutate(term = factor(term, levels = unique(ROC_long$term))) %>%
  arrange(term) %>% # easier plotting later
  left_join(ROC_tbl_scores %>% select(term, shortname), by='term')
ROC_youden_RS <- ROC_long_RS %>%
  left_join(ROC_youden, by='term') %>%
  filter(youden_J == max_youden_J)

# plot
gg_mm2 <- ggplot(ROC_long_RS, aes(x = 1 - specificity, y = sensitivity,
                              color = shortname)) +
  geom_line(aes(group=term)) +
  geom_point(data=ROC_youden_RS) +
  geom_label_repel(data=ROC_youden_RS, aes(label=shortname), size=2)
gg_mm2 <- add_ROC_plot_elements(gg, legend=FALSE)
gg_mm2 <- gg_mm2 + theme(
  axis.title = element_text(size=7),
  axis.text = element_text(size=6))
loc_fig <- paste0(dir_figs,"main_models_ROCs")
ggsave(paste0(loc_fig,".png"), gg_mm2, width=180, height=180, units="mm", dpi=300)


ggarrange(plotlist = list(gg_mm2, gg_mm1), ncol=2, labels='AUTO',
          align='h', widths=c(1.5,1), vjust=2)
loc_fig <- paste0(dir_figs,"main_models_combined")
ggsave(paste0(loc_fig,".png"), width=240, height=140, units="mm", dpi=300)

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
