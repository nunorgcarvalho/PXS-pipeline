library(tidyverse)
library(data.table)
source('code/00_paths.R')

pheno <- as_tibble(fread('scratch/cohorts/BRS_cohort_ALL.txt'))

pheno %>% select(overweight, `BRS-normal_BMI-cov_bvr`, `BRS-high_BMI-cov_bvr`) %>%
  drop_na() %>%
  group_by(overweight) %>%
  summarize(N=n(),
            `BRS-normal_BMI-cov_bvr.MEAN` = mean(`BRS-normal_BMI-cov_bvr`),
            `BRS-normal_BMI-cov_bvr.SE` = sd(`BRS-normal_BMI-cov_bvr`) / sqrt(N),
            `BRS-high_BMI-cov_bvr.MEAN` = mean(`BRS-high_BMI-cov_bvr`),
            `BRS-high_BMI-cov_bvr.SE` = sd(`BRS-high_BMI-cov_bvr`) / sqrt(N) )

## BMI comparisons ####
loc_pheno_full1 <- "/n/no_backup2/patel/uk_biobank/main_data_34521/ukb34521.tab"
cols_to_keep <- c('f.eid', paste0('f.21001.',0:3,'.0') )
all_BMIs <- as_tibble(fread(loc_pheno_full1, select=cols_to_keep))

pheno2 <- pheno %>% left_join(
  all_BMIs %>%
    rename(IID=f.eid, BMI_t0 = `f.21001.0.0`,
           BMI_t1 = `f.21001.1.0`, BMI_t2 = `f.21001.2.0`), by='IID') %>%
  select(IID, `BRS-ALL-cov_bvr`, `CRS-ALL`, `PRS_T2D`,
         starts_with('T2D_onset'), overweight, starts_with('BMI_t')) %>%
  # pivot_longer(cols=starts_with('BMI_t'), names_to='time',
  #              names_prefix='BMI_t', values_to='BMI') %>%
  mutate(T2D_onset=as.factor(T2D_onset),
         overweight=as.factor(overweight))

pheno3 <- pheno2 %>%
  mutate(BMI_t0.1 = BMI_t1 - BMI_t0,
         BMI_t0.2 = BMI_t2 - BMI_t0,
         BMI_change = ifelse(!is.na(BMI_t0.2), BMI_t0.2, BMI_t0.1))

dBMI_tbl <- tibble(expand.grid(T2D=c(0,1),overweight=c(0,1))) %>%
  mutate(dBMI_mean=NA, dBMI_sd=NA, N=NA,
         BRS_cor = NA, CRS_cor = NA, PRS_cor = NA,
         BRS_cor_se = NA, CRS_cor_se = NA, PRS_cor_se = NA,
         BRS_cor_p = NA, CRS_cor_p = NA, PRS_cor_p = NA )

get_se <- function(cor, CI, level=0.95) {
  CI_z <- qnorm(1 - ((1-level) / 2))
  cor_se <- (cor1$conf.int[2] - cor1$conf.int[1]) / (CI_z*2)
  return(cor_se)
}

col_risk_score <- list('BRS' = 'BRS-ALL-cov_bvr',
                       'CRS' = 'CRS-ALL',
                       'PRS' = 'PRS_T2D')

for (i in 1:nrow(dBMI_tbl)) {
  pheno_subset <- pheno3 %>%
    filter(T2D_onset == dBMI_tbl$T2D[i],
           overweight == dBMI_tbl$overweight[i]) %>%
    filter(!is.na(BMI_change))
  
  for (risk_score in c('BRS','CRS',"PRS")) {
    cor1 <- cor.test(pheno_subset$BMI_change,
                     pheno_subset[[ col_risk_score[[risk_score]] ]])
    dBMI_tbl[i,paste0(risk_score,'_cor') ] <- cor1$estimate
    dBMI_tbl[i,paste0(risk_score,'_cor_se') ] <- get_se(cor1$conf.int)
    dBMI_tbl[i,paste0(risk_score,'_cor_p') ] <- cor1$p.value
  }
  
  dBMI_tbl$dBMI_mean[i] <- mean(pheno_subset$BMI_change, na.rm=TRUE)
  dBMI_tbl$dBMI_sd[i] <- sd(pheno_subset$BMI_change, na.rm=TRUE)
  dBMI_tbl$N[i] <- nrow(pheno_subset)
  
}

CI95_z <- qnorm(1 - ((1-0.95) / (2*4)))
# long code, im sorry
dBMI_tbl_long <- dBMI_tbl %>%
  pivot_longer(cols = ends_with('_cor'),
               names_to = "risk_score", 
               values_to = "cor") %>%
  mutate(risk_score = str_remove(risk_score, "_cor")) %>%
  left_join(
    dBMI_tbl %>%
      pivot_longer(cols = ends_with('_cor_se'), 
                   names_to = "risk_score", 
                   values_to = "cor_se") %>%
      mutate(risk_score = str_remove(risk_score, "_cor_se")),
    by = c("T2D", "overweight", "risk_score")
  ) %>%
  
  left_join(
    dBMI_tbl %>%
      pivot_longer(cols = ends_with('_cor_p'), 
                   names_to = "risk_score", 
                   values_to = "cor_p") %>%
      mutate(risk_score = str_remove(risk_score, "_cor_p")),
    by = c("T2D", "overweight", "risk_score")
  ) %>%
  select(T2D, overweight, risk_score, cor, cor_se, cor_p) %>%
  mutate(cor_low = cor - CI95_z*cor_se,
         cor_upp = cor + CI95_z*cor_se,
         T2D_label = ifelse(T2D==1,'T2D','no T2D'),
         OW_label = ifelse(overweight==1,'high BMI','normal BMI'))




ggplot(dBMI_tbl_long, aes(x=risk_score, y=cor)) +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_point() +
  geom_errorbar(aes(ymin=cor_low, ymax=cor_upp), width=0.25) +
  facet_grid(T2D_label ~ OW_label) +
  labs(x = 'Risk Score',
       y = 'Correlation with change in BMI over time',
       subtitle = 'Adjusted for multiple testing (4 groups)') +
  theme_bw()
