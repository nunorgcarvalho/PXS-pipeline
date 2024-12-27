# computes coefficients for T2D-BRS, assesses accuracy, and updates pheno file

## Libraries and directories ##
library(tidyverse)
library(data.table)
source('code/paths.R')

fields <- as_tibble(fread('scratch/fields_tbl.txt'))
col_covs <- fields$term[fields$use_type == 'covar']
col_behaviors <- fields$term[fields$use_type == 'behavior']

pheno <- as_tibble(fread('scratch/pheno.txt'))

# cross validation cox LASSO regression ####
get_cvglm_obj <- function(training_tbl, cols_use) {
  glm_y <- survival::Surv(training_tbl$T2D_onset_days,training_tbl$T2D_onset)
  glm_x <- as.matrix( training_tbl[,cols_use] )
  
  cv_glm1 <- glmnet::cv.glmnet(glm_x, glm_y, family='cox', alpha=1, nfolds=10,
                               type.measure = 'C')
  return(cv_glm1)
}

# BMI cohorts ####
cohorts <- list(
  'ALL' = pheno,
  'normal_BMI' = pheno %>% filter(overweight==0),
  'high_BMI' = pheno %>% filter(overweight==1)
)
lasso_factors <- list(
  'cov' = c(col_covs),
  'cov_bvr' = c(col_covs,col_behaviors),
  'cov_BMI' = c(col_covs, 'f21001'),
  'cov_BMI_bvr' = c(col_covs, 'f21001', col_behaviors)
)

set.seed(2024)
models <- list()
for (i in 1:length(cohorts)) {
  cohort <- names(cohorts)[i]
  training_tbl <- cohorts[[i]] %>% filter(sample_group=='A')
  
  for (j in 1:length(lasso_factors)) {
    training_tbl_j <- training_tbl %>% 
      select(T2D_onset, T2D_onset_days,
             all_of(lasso_factors[[j]])) %>%
      drop_na()
    factor_set <- names(lasso_factors)[j]
    
    model_name <- paste0(cohort,'-',factor_set)
    print(paste(i, j, model_name))
    models[[model_name]] <- get_cvglm_obj(training_tbl_j, lasso_factors[[j]])
    print(paste(i, j, model_name,models[[model_name]]$cvm[models[[1]]$index[2]]))
  }
}

dir.create('scratch/BRS_models', showWarnings = FALSE)
save(models, file='scratch/BRS_models/cv_glm_models.RData')
#load('scratch/BRS_models/cv_glm_models.RData')

# computes each BRS for all individuals ####
# also saves coefficients
# and also computes C_statisc
BRS_coeffs <- tibble(term='')[0,]
tbl_Cstat <- tibble(training_cohort='', model_factors='',testing_cohort='',
                    model_name='',model_label='', C_stat=0,C_stat_se=0)[0,]
for (i in 1:length(models)) {
  cv_glm <- models[[i]]
  col_BRS <- paste0('BRS-',names(models)[i])
  
  coef_matrix <- as.matrix(cv_glm$glmnet.fit$beta)
  index_1se <- cv_glm$index[2] #1se
  betas <- coef_matrix[,index_1se]
  
  print(paste(i,'saving betas'))
  
  col_beta <- paste0('beta-',col_BRS)
  if (!exists('BRS_coeffs')) {
    BRS_coeffs$term <- names(betas)
    BRS_coeffs[[col_beta]] <- betas
  } else {
    betas_tbl <- as_tibble(betas, rownames='term')
    colnames(betas_tbl)[2] <- col_beta
    BRS_coeffs <- BRS_coeffs %>%
      full_join(betas_tbl)
  }
  
  
  ## computing and C-statistic ####
  model_name <- names(models)[i]
  print(model_name)
  model_settings <- str_split(model_name,'-')[[1]]
  new_col_BRS <- paste0('BRS-',model_settings[2])
  
  col_BRS <- paste0('BRS-',model_name)
  col_beta <- paste0('beta-',col_BRS)
  formula <- as.formula(paste0(
    'survival::Surv(T2D_onset_days, T2D_onset) ~ `',col_BRS,'`') )
  
  
  for (j in 1:length(cohorts)) {
    testing_group_name <- names(cohorts)[j]
    print(paste(i,j))
    factors_matrix <- as.matrix( cohorts[[j]][,names(betas)] )
    BRS_raw <- (factors_matrix %*% betas)[,1]
    BRS_vec <- scale(BRS_raw)[,1]
    
    cohorts[[j]][[col_BRS]] <- BRS_vec
    
    if (model_settings[1] == testing_group_name & j != 1) {
      if ( !(new_col_BRS %in% colnames(cohorts[['ALL']])) ) {
        cohorts[['ALL']][[new_col_BRS]] <- NA_real_
      }
      
      cohort_indices <- cohorts[['ALL']]$IID %in% cohorts[[j]]$IID
      cohorts[['ALL']][cohort_indices,new_col_BRS] <- BRS_vec
    }
    
    data_testing <- cohorts[[j]] %>% filter(sample_group=='B')
    sc1 <- survival::concordance(formula, data=data_testing, reverse=TRUE)
    
    tbl_Cstat <- tbl_Cstat %>% add_row(
      training_cohort = model_settings[1],
      model_factors = model_settings[2],
      testing_cohort = testing_group_name, 
      model_name = model_name,
      model_label = str_replace_all(model_factors, '_',' + '),
      C_stat = sc1$concordance, C_stat_se = sqrt(sc1$var) )
    
  }
}
BRS_coeff_table <- fields %>%
  select(term,traitname) %>%
  filter(term %in% BRS_coeffs$term) %>%
  left_join(BRS_coeffs)
fwrite(BRS_coeff_table, 'scratch/BRS_models/BRS_coefficients.txt', sep='\t')
fwrite(tbl_Cstat, 'scratch/BRS_models/BRS_Cstats.txt')

# saves cohorts w/ BRS to system
dir.create('scratch/cohorts/', showWarnings = FALSE)
for (j in 1:length(cohorts)) {
  cohort_name <- names(cohorts)[j]
  cohort_tbl <- cohorts[[j]]
  
  loc_out <- paste0('scratch/cohorts/BRS_cohort_',cohort_name,'.txt')
  fwrite(cohort_tbl, file=loc_out, sep='\t')
}

## visualizing C-stat ####
ci95_z <- -1* qnorm((1 - 0.95)/2)
ggplot(tbl_Cstat, aes(x = testing_cohort , y=C_stat,
                      color = training_cohort)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = C_stat - ci95_z * C_stat_se,
                    ymax = C_stat + ci95_z * C_stat_se),
                width= 0.5, position = position_dodge(width = 0.5)) +
  scale_y_continuous(breaks = seq(0,1,1/100),
                     minor_breaks = seq(0,1,1/400)) +
  theme_bw() +
  facet_wrap(~model_label, nrow=1) +
  theme(legend.position = 'top')

