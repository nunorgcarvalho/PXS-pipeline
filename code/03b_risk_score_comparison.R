# compares C-statistics for BRS, CRS, and PRS

## Libraries and directories ##
library(tidyverse)
library(data.table)
source('code/00_paths.R')

# reads cohort data
cohort_tbl <- as_tibble(fread(paste0(dir_scratch,'cohorts/BRS_cohort_ALL.txt')))

scores <- list(
  'cov' = c('BRS-ALL-cov'),
  'BRS' = c('BRS-ALL-cov_bvr'),
  'CRS' = c('CRS-ALL'),
  'PRS' = c('PRS_T2D'),
  'BRS+CRS' = c('BRS-ALL-cov_bvr','CRS-ALL_nocov'),
  'BRS+PRS' = c('BRS-ALL-cov_bvr','PRS_T2D'),
  'CRS+PRS' = c('CRS-ALL','PRS_T2D'),
  'BRS+CRS+PRS' = c('BRS-ALL-cov_bvr','CRS-ALL_nocov','PRS_T2D')
)

# combines CRS, PRS, BRS using linear weights

model_performances <- tibble(model_name='',C_stat=0, C_stat_se=0)[0,]

for (i in 1:length(scores)) {
  model_name = names(scores)[i]
  model_cols = scores[[i]]
  
  print(paste(i,model_name))
  
  if (length(model_cols) > 1) {
    # determines linear weight of each factor
    glm_tbl <- cohort_tbl %>%
      filter(sample_group == 'A') %>%
      select(T2D_onset_days, T2D_onset, all_of(model_cols) ) %>%
      drop_na()
    glm_y <- survival::Surv(glm_tbl$T2D_onset_days,glm_tbl$T2D_onset)
    glm_x <- as.matrix( glm_tbl[,model_cols] )
    
    formula <- as.formula(paste0(
      'survival::Surv(T2D_onset_days, T2D_onset) ~ .') )
    glm1 <- glmnet::glmnet(glm_x, glm_y, family='cox', lambda=0, data=glm_tbl)
    betas <- as.matrix(glm1$beta)[,1]
    
    factors_matrix <- as.matrix( cohort_tbl[,names(betas)] )
    score_raw <- (factors_matrix %*% betas)[,1]
    score_vec <- scale(score_raw)[,1]
  } else {
    score_vec <- cohort_tbl[[model_cols[1]]]
  }
  cohort_tbl[[model_name]] <- score_vec
  
  formula <- as.formula(paste0(
    'survival::Surv(T2D_onset_days, T2D_onset) ~ `',model_name,'`') )
  sc1 <- survival::concordance(formula, reverse=TRUE,
                               data=cohort_tbl %>% filter(sample_group=='B'))
  
  model_performances <- model_performances %>% add_row(
    model_name = model_name,
    C_stat = sc1$concordance, C_stat_se = sqrt(sc1$var)
  )
  
}

model_performances$model_name <- factor(model_performances$model_name,
                                        levels = model_performances$model_name)
# plots performance

ggplot(model_performances, aes(x=fct_rev(model_name))) +
  geom_point(aes(y=C_stat)) +
  geom_errorbar(aes(ymin = C_stat - 1.96*C_stat_se,
                    ymax = C_stat + 1.96*C_stat_se),
                width=0.5) +
  theme_bw() +
  labs(x = 'Model Components',
       y = 'C-statistic (95% CI)') +
  coord_flip()
