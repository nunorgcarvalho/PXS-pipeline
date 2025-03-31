# computes coefficients for T2D-BRS, assesses accuracy, and updates pheno file

## Libraries and directories ##
library(tidyverse)
library(data.table)
source('code/00_paths.R')

fields <- as_tibble(fread(paste0(dir_scratch,'fields_tbl.txt')))
col_covs <- fields$term[fields$use_type == 'covar']
col_behaviors <- fields$term[fields$use_type == 'behavior']
col_SES <- fields$term[fields$use_type == 'SES']

pheno <- as_tibble(fread(paste0(dir_scratch,'pheno.txt')))

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
  'high_BMI' = pheno %>% filter(overweight==1),
  'females' = pheno %>% filter(sex==0),
  'males' = pheno %>% filter(sex == 1),
  'T2D_early' = pheno %>% filter(T2D_onset_cat %in% c('early','control')),
  'T2D_late' = pheno %>%  filter(T2D_onset_cat %in% c('late','control'))
)
lasso_factors <- list(
  'cov' = c(col_covs),
  'cov_bvr' = c(col_covs,col_behaviors),
  'cov_BMI' = c(col_covs, 'f21001'),
  'cov_BMI_bvr' = c(col_covs, 'f21001', col_behaviors),
  'cov_SES' = c(col_covs,col_SES),
  'cov_bvr_SES' = c(col_covs,col_behaviors,col_SES)
)

# makes table of models to run LASSO on
models2run <- expand_grid(cohort = names(cohorts)[1:3], factor_set = names(lasso_factors)[1:4]) %>%
  as_tibble() %>%
  add_row(expand_grid(cohort = names(cohorts)[1], factor_set = names(lasso_factors)[5:6])) %>%
  add_row(expand_grid(cohort = names(cohorts)[4:7], factor_set = names(lasso_factors)[1:2] ))

set.seed(2024)
models <- list()
######
for (i in 1:nrow(models2run)) {
  cohort <- models2run$cohort[i]
  factor_set <- models2run$factor_set[i]
  model_name <- paste0(cohort,'-',factor_set)
  if (model_name %in% names(models)) {
    print(paste(i, model_name,'SKIPPING'))
    next
  }
  print(paste(i, model_name))
  
  training_tbl_j <- cohorts[[cohort]] %>% filter(sample_group=='A') %>%
    select(T2D_onset, T2D_onset_days,
           all_of(lasso_factors[[factor_set]])) %>%
    drop_na()
  
  models[[model_name]] <- get_cvglm_obj(training_tbl_j, lasso_factors[[factor_set]])
  print(paste(i, model_name, models[[model_name]]$cvm[models[[model_name]]$index[2]]))
  
  # gets standard deviations of factors, for later rescaling if desired
  models[[model_name]]$factor_SDs <- sapply(training_tbl_j[,-c(1,2)], sd)
  
  # makes BMI-adjusted cov-bvr models
  if (factor_set == 'cov_BMI_bvr') {
    # BMIadj1: removes direct effect of BMI on T2D
    # by setting model coefficient to 0
    full_model <- models[[model_name]]
    BMI_index <- which(full_model$glmnet.fit$beta@Dimnames[[1]] == 'f21001')
    full_model$glmnet.fit$beta[BMI_index,] <- 0
    
    model_name <- paste0(cohort,'-',str_replace(factor_set,'_BMI',''), '_BMIadj1')
    print(paste(i, model_name))
    models[[model_name]] <- full_model
  }
  
}

dir.create(paste0(dir_scratch,'BRS_models'), showWarnings = FALSE)
#save(models, file=paste0(dir_scratch,'BRS_models/cv_glm_models.RData'))
#load(paste0(dir_scratch,'BRS_models/cv_glm_models.RData'))

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
    print(paste(i,j, testing_group_name))
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
    
    # makes BMIadj version of cov_bvr model
    if (model_settings[2] == 'cov_BMI_bvr') {
      # BMIadj2: removes all association of BMI on T2D
      # by regressing the cov_BMI_bvr score on BMI, and taking the resid
      BMIadj2_formula <- as.formula(paste0('`',col_BRS, '` ~ f21001') )
      # BMI-adjustment learned in testing cohort to guarantee independence during
      # C-statistic calculation later
      lm1_BRS_BMI <- lm(BMIadj2_formula, data=data_testing) 
      BRS_BMI_pred <- predict(lm1_BRS_BMI, newdata=cohorts[[j]])
      adj2_BRS_raw <- cohorts[[j]][[col_BRS]] - BRS_BMI_pred
      adj2_BRS_vec <- scale(adj2_BRS_raw)[,1]
      
      adj2_col_BRS <- paste0(str_replace(col_BRS,'_BMI',''), '_BMIadj2')
      cohorts[[j]][[adj2_col_BRS]] <- adj2_BRS_vec
      
      data_testing <- cohorts[[j]] %>% filter(sample_group=='B')
      adj2_formula <- as.formula(paste0(
        'survival::Surv(T2D_onset_days, T2D_onset) ~ `',adj2_col_BRS,'`') )
      sc2 <- survival::concordance(adj2_formula, data=data_testing, reverse=TRUE)
      
      tbl_Cstat <- tbl_Cstat %>% add_row(
        training_cohort = model_settings[1],
        model_factors = paste0(str_replace(model_settings[2],'_BMI',''), '_BMIadj2'),
        testing_cohort = testing_group_name, 
        model_name = paste0(str_replace(model_name,'cov_BMI','cov_'), '_BMIadj2'),
        model_label = str_replace_all(model_factors, '_',' + '),
        C_stat = sc2$concordance, C_stat_se = sqrt(sc2$var) )
    }
    
  }
}
tbl_Cstat$model_label <- str_replace(tbl_Cstat$model_label,'\\+ BMIadj1','(BMI-adj1)')
tbl_Cstat$model_label <- str_replace(tbl_Cstat$model_label,'\\+ BMIadj2','(BMI-adj2)')

BRS_coeff_table <- fields %>%
  select(term,traitname) %>%
  filter(term %in% BRS_coeffs$term) %>%
  left_join(BRS_coeffs)
fwrite(BRS_coeff_table, 'scratch/BRS_models/BRS_coefficients.txt', sep='\t')
fwrite(tbl_Cstat, 'scratch/BRS_models/BRS_Cstats.txt')

# saves cohorts w/ BRS to system
dir.create(paste0(dir_scratch,'cohorts/'), showWarnings = FALSE)
for (j in 1:length(cohorts)) {
  cohort_name <- names(cohorts)[j]
  cohort_tbl <- cohorts[[j]]
  
  loc_out <- paste0(dir_scratch,'cohorts/BRS_cohort_',cohort_name,'.txt')
  fwrite(cohort_tbl, file=loc_out, sep='\t', na="NA", quote=FALSE)
}

## visualizing C-stat ####
ci95_z <- -1* qnorm((1 - 0.95)/2)
tbl_Cstat <- tbl_Cstat %>%
  mutate(C_stat_low = C_stat - ci95_z * C_stat_se,
         C_stat_upp = C_stat + ci95_z * C_stat_se)
ggplot(tbl_Cstat, aes(x = testing_cohort , y=C_stat,
                      color = training_cohort)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = C_stat_low,
                    ymax = C_stat_upp),
                width= 0.5, position = position_dodge(width = 0.5)) +
  scale_y_continuous(breaks = seq(0,1,1/100),
                     minor_breaks = seq(0,1,1/400)) +
  theme_bw() +
  facet_wrap(~model_label, nrow=1) +
  theme(legend.position = 'top')
