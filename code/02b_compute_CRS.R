# computes clinical risk score (CRS) for T2D using 6 clinically relevant factors
# plus covariates. Very similar code to 02a_compute_BRS.R

## Libraries and directories ##
library(tidyverse)
library(data.table)
source('code/00_paths.R')

fields <- as_tibble(fread('scratch/fields_tbl.txt'))
col_covs <- fields$term[fields$use_type == 'covar']
col_CRFs <- fields$term[fields$use_type == 'CRF']

pheno <- as_tibble(fread('scratch/pheno.txt')) #%>%
  # rename(
  #   'SBP' = 'f4080',
  #   'BMI' = 'f21001',
  #   'Glucose' = 'f30740',
  #   'HbA1c' = 'f30750',
  #   'HDL' = 'f30760',
  #   'Triglycerides' = 'f30870'
  # )

# loads list of CRFs
loc_CFRs_tbl <- paste0(dir_script,"../input_data/CRFs_table.txt")
CRFs_tbl <- as_tibble(fread(loc_CFRs_tbl))

# cross validation cox LASSO regression ####
get_cvglm_obj <- function(training_tbl, cols_use) {
  glm_y <- survival::Surv(training_tbl$T2D_onset_days,training_tbl$T2D_onset)
  glm_x <- as.matrix( training_tbl[,cols_use] )
  
  cv_glm1 <- glmnet::cv.glmnet(glm_x, glm_y, family='cox', alpha=1, nfolds=10,
                               type.measure = 'C')
  return(cv_glm1)
}

cols_use <- c(col_covs, col_CRFs[-4])
training_tbl <- pheno %>% filter(sample_group == 'A') %>%
  select(T2D_onset, T2D_onset_days, all_of(cols_use)) %>%
  drop_na()

cvglm_CRS <- get_cvglm_obj(training_tbl, cols_use)

coef_matrix <- as.matrix(cvglm_CRS$glmnet.fit$beta)
index_1se <- cvglm_CRS$index[2] #1se
betas <- coef_matrix[,index_1se]

compute_score <- function(betas) {
  factors_matrix <- as.matrix( cohort_tbl[,names(betas)] )
  score_raw <- (factors_matrix %*% betas)[,1]
  score_vec <- scale(score_raw)[,1]
  return(score_vec)
}

# read cohorts table
loc_cohort <- paste0(dir_scratch,'cohorts/BRS_cohort_ALL.txt')
cohort_tbl <- as_tibble(fread(loc_cohort))

# computes CRS and saves it, excludes covariate scores
cohort_tbl$`CRS-ALL_nocov` <- compute_score(betas)
betas[names(betas) %in% col_covs] <- 0
cohort_tbl$`CRS-ALL_nocov` <- compute_score(betas)

fwrite(cohort_tbl, file=loc_cohort, sep='\t', na="NA", quote=FALSE)

# # computes C-statistic
# data_testing <- cohort_tbl %>% filter(sample_group=='B') %>%
#   mutate(`CRS+BRS` = `CRS-ALL_nocov` + `BRS-ALL-cov_bvr`)
# formula <- as.formula(paste0(
#   'survival::Surv(T2D_onset_days, T2D_onset) ~ `CRS+BRS`') )
# sc1 <- survival::concordance(formula, data=data_testing, reverse=TRUE)
# sc1
