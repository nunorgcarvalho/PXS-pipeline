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
cols_to_keep <- c('f.eid', paste0('f.21001.',0:2,'.0') )
all_BMIs <- as_tibble(fread(loc_pheno_full1, select=cols_to_keep))

pheno2 <- pheno %>% left_join(
  all_BMIs %>% rename(IID=f.eid, BMI_t0 = `f.21001.0.0`,
                      BMI_t1 = `f.21001.1.0`, BMI_t2 = `f.21001.2.0`,),
  by='IID'
) %>% select(IID,starts_with('T2D_onset'), overweight, starts_with('BMI_t')) %>%
  filter(!is.na(BMI_t2)) %>%
  pivot_longer(cols=starts_with('BMI_t'), names_to='time',
               names_prefix='BMI_t', values_to='BMI') %>%
  mutate(T2D_onset=as.factor(T2D_onset), overweight=as.factor(overweight)) %>%
  filter(time != 1, !is.na(overweight))
  

ggplot(pheno2, aes(x=time, y=BMI, color=T2D_onset,
                   linetype=overweight)) +
  geom_boxplot(alpha=0.1, shape=1) +
  #facet_wrap(~ overweight) +
  theme_bw()
