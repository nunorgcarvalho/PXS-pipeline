# Looks at how BRS changes with sex-specific effects, time-specific effects, and SES

# Libraries and directories ####
library(tidyverse)
library(data.table)
library(pROC)
library(ggpubr)
source('code/00_paths.R')
source('code/00_plotting.R')

pheno <- as_tibble(fread(paste0(dir_scratch,'cohorts/BRS_cohort_ALL.txt'))) %>%
  mutate(SES = factor(f738, levels = c(1,2,3,4,5),
               labels = c('1' = '<£18k',
                          '2' = '£18k-£31k',
                          '3' = '£31k-£52k',
                          '4' = '£52k-£100k',
                          '5' = '>£100k')))
test <- pheno %>% filter(sample_group == 'B') %>%
  mutate(sex = factor(sex, levels=c(0,1),
                      labels = c('Female', 'Male')))
         
# shared scatterplot theme stuff ####
add_sensitivity_elements <- function(gg, cor_obj, legend_top=TRUE) {
  gg <- gg +
    geom_point(alpha=0.1, shape=1) +
    geom_abline(slope=1, color=dash_color, linetype='dashed') +
    geom_smooth(method='lm', formula='y ~ x', se=FALSE, size=0.5) +
    annotate('text', label = annotate_cor.test(cor_obj), # throws out harmless warning message
             x=Inf, y=-Inf, vjust=-0.5, hjust=1.05, parse=FALSE,
             size=3) +
    coord_fixed(expand = FALSE) +
    theme(axis.title = element_text(size=8),
          axis.text = element_text(size=6),
          legend.title = element_text(size=6),
          legend.text = element_text(size=5))
  if (legend_top) {gg <- gg + theme(
    legend.position = 'top',
    legend.margin = margin(b = -10) )}
  
  return(gg)
}

# sex-specific analysis ####
cor_sex <- cor.test(test$`BRS-females-cov_bvr`,test$`BRS-males-cov_bvr`)
gg_sex <- ggplot(test, aes(x=`BRS-females-cov_bvr`,y=`BRS-males-cov_bvr`,
                  color = sex)) +
  labs(x = 'BRS (trained in females)',
       y = 'BRS (trained in males)',
       color = 'Sex')
gg_sex <- add_sensitivity_elements(gg_sex, cor_sex)
#loc_fig <- paste0(dir_figs,"BRS_sex_specific_comparison")
#ggsave(paste0(loc_fig,".png"), width=180, height=180, units="mm", dpi=300)

# time-specific analysis ####
cor_time <- cor.test(test$`BRS-T2D_early-cov_bvr`,test$`BRS-T2D_late-cov_bvr`)
gg_time <- ggplot(test, aes(x=`BRS-T2D_early-cov_bvr`,y=`BRS-T2D_late-cov_bvr`,
                 color = T2D_onset_cat)) +
  labs(x = 'BRS (trained with early T2D diagnoses)',
       y = 'BRS (trained with late T2D diagnoses)',
       color = 'T2D Diagnosis')
gg_time <- add_sensitivity_elements(gg_time, cor_time)
#loc_fig <- paste0(dir_figs,"BRS_time_specific_comparison")
#ggsave(paste0(loc_fig,".png"), width=180, height=180, units="mm", dpi=300)

## combined plot ####
ggarrange(plotlist = list(gg_sex, gg_time), ncol=2, labels='AUTO',
          align='h', widths=c(1,1))
loc_fig <- paste0(dir_figs,"BRS_sensitivity_sex_time")
ggsave(paste0(loc_fig,".png"), width=180, height=100, units="mm", dpi=300)

# SES analysis ####

## coefficient comparison ####
load(paste0(dir_scratch,'BRS_models/cv_glm_models.RData'))
model_bvr_SES <- models[['ALL-cov_bvr_SES']]
bSES_bvr_SES <- coef(model_bvr_SES, s='lambda.1se')[length(model_bvr_SES$factor_SDs)] *
  model_bvr_SES$factor_SDs[length(model_bvr_SES$factor_SDs)]
model_SES <- models[['ALL-cov_SES']]
bSES_SES <- coef(model_SES, s='lambda.1se')[length(model_SES$factor_SDs)] *
  model_SES$factor_SDs[length(model_SES$factor_SDs)]
bSES_bvr_SES %>% exp()
bSES_SES %>% exp()

## regression coeffs ####
# runs regression on BRS (w/ SES) ~ BRS (w/o SES) for each SES group
SES_lm <- tibble(SES = levels(test$SES), a = NA,a_se = NA, b = NA, b_se=NA)
for (i in 1:length(levels(test$SES))) {
  SES_group <- levels(test$SES)[i]
  test_SES <- test %>% filter(SES == SES_group)
  lm_SES <- lm(`BRS-ALL-cov_bvr_SES` ~ `BRS-ALL-cov_bvr`, data=test_SES )
  SES_lm$a[i] <- summary(lm_SES)$coefficients[1,1]
  SES_lm$a_se[i] <- summary(lm_SES)$coefficients[1,2]
  SES_lm$b[i] <- summary(lm_SES)$coefficients[2,1]
  SES_lm$b_se[i] <- summary(lm_SES)$coefficients[2,2]
}
SES_lm <- SES_lm %>% mutate(
  a_low = a - CI95_z * a_se,
  a_upp = a + CI95_z * a_se,
  b_low = b - CI95_z * b_se,
  b_upp = b + CI95_z * b_se
)
gg_SES_lm <- ggplot(SES_lm, aes(x=a,y=b, color=SES)) +
  geom_point() +
  geom_text(aes(label = SES, y=b_upp+0.002)) +
  geom_vline(xintercept=0, color=dash_color, linetype='dashed') +
  geom_hline(yintercept=1, color=dash_color, linetype='dashed') +
  geom_errorbar( aes(ymin=b_low, ymax=b_upp), width =0.000) +
  geom_errorbarh(aes(xmin=a_low, xmax=a_upp), height=0.000) +
  scale_x_continuous(expand = c(0.1,0)) +
  labs(x = 'Intercept term of (BRS w/ SES) ~ (BRS w/o SES)',
       y = 'Beta term of (BRS w/ SES) ~ (BRS w/o SES)',
       color = 'Income') +
  theme(legend.position='none',
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))
#loc_fig <- paste0(dir_figs,"BRS_SES_lm")
#ggsave(paste0(loc_fig,".png"), width=150, height=180, units="mm", dpi=300)

## testing interaction ####
SES_models <- tibble(
  factor_set = c('cov','cov_bvr','cov_bvr_SES','cov_SES','BRS * SES','BRS (w/ SES) * SES'),
  AUC = NA, AUC_low = NA, AUC_upp = NA
)
SES_ROCs <- list()

pheno_SES <- pheno %>%
  select(T2D_onset_days, T2D_onset, sample_group, f738,
         any_of(paste0('BRS-ALL-',SES_models$factor_set))) %>%
  drop_na()
colnames(pheno_SES) <- str_replace_all(colnames(pheno_SES), 'BRS-ALL-','')
train_SES <- pheno_SES %>% filter(sample_group == 'A')
test_SES <- pheno_SES %>% filter(sample_group == 'B')
test_y <- test_SES$T2D_onset

for (i in 1:nrow(SES_models)) {
  factor_set <- SES_models$factor_set[i]
  print(paste(i, factor_set))
  if (factor_set == 'BRS * SES') {
    sv1 <- survival::coxph(survival::Surv(T2D_onset_days, T2D_onset) ~ `cov_bvr` * f738,data=train_SES)
    test_SES[[factor_set]] <- ( predict(sv1, newdata=test_SES) %>% scale() )[,1]
  } else if (factor_set == 'BRS (w/ SES) * SES') {
    sv1 <- survival::coxph(survival::Surv(T2D_onset_days, T2D_onset) ~ `cov_bvr_SES` * f738,data=train_SES)
    test_SES[[factor_set]] <- ( predict(sv1, newdata=test_SES) %>% scale() )[,1]
  } 
  pred_vec <- test_SES[[factor_set]]
  roc1 <- roc(test_y, pred_vec, ci=TRUE)
  SES_ROCs[[factor_set]] <- roc1
}
SES_models <- SES_models %>% mutate(
  AUC = map_dbl(SES_ROCs, 'auc'), # extracts AUC
  AUC_low = map_dbl(SES_ROCs, ~ as.numeric(.$ci)[1] ),
  AUC_upp = map_dbl(SES_ROCs, ~ as.numeric(.$ci)[3] ),
  #model_label = str_replace_all(factor_set, '_',' + ')
  model_label = factor(
    factor_set,
    levels = c('cov','cov_SES','cov_bvr','cov_bvr_SES', 'BRS * SES', 'BRS (w/ SES) * SES'),
    labels = c('Covariates-only', 'Covariates + SES', 'Behavioral Risk Score (BRS)',
               'BRS + SES', 'BRS * SES', '(BRS + SES) * SES'))
) %>% arrange(model_label)

gg_SES_AUC <- ggplot(SES_models, aes(x=fct_rev(model_label), y=AUC)) +
  geom_point() +
  geom_errorbar(aes(ymin=AUC_low, ymax=AUC_upp), width=0.25) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0,1,1/20), minor_breaks = seq(0,1,1/100)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width=15)) +
  labs(x = 'Model',
       y = 'AUC in testing cohort with SES data') +
  theme(legend.position='none',
        axis.title = element_text(size=7),
        axis.text = element_text(size=6))
#loc_fig <- paste0(dir_figs,"BRS_SES_AUCs")
#ggsave(paste0(loc_fig,".png"), width=150, height=180, units="mm", dpi=300)

### combined plots ####
ggarrange(plotlist = list(gg_SES_lm, gg_SES_AUC), ncol=2, labels='AUTO',
          align='h', widths=c(1.5,1), vjust=2)
loc_fig <- paste0(dir_figs,"BRS_SES2")
ggsave(paste0(loc_fig,".png"), width=240, height=140, units="mm", dpi=300)


## DeLong test for comparing ROCs ####
roc.test(SES_ROCs[['cov_bvr']], SES_ROCs[['cov_SES']])
roc.test(SES_ROCs[['cov_bvr']], SES_ROCs[['cov_bvr_SES']])
roc.test(SES_ROCs[['cov_bvr']], SES_ROCs[['BRS * SES']])
roc.test(SES_ROCs[['cov_bvr']], SES_ROCs[['BRS (w/ SES) * SES']])


## BRS comparison ####
cor_SES <- cor.test(test$`BRS-ALL-cov_bvr`,test$`BRS-ALL-cov_bvr_SES`)
cor_SES2 <- cor.test(test$`BRS-ALL-cov_bvr`,test$`BRS-ALL-cov_SES`)
gg_SES <- ggplot(test, aes(x=`BRS-ALL-cov_bvr`,y=`BRS-ALL-cov_bvr_SES`,
                           color = SES)) +
  labs(x = 'BRS (without SES factor)',
       y = 'BRS (with SES factor)',
       color = 'Income') +
  theme(legend.position = c(0,1), legend.justification = c(-0.05, 1.05),
        legend.box.background = element_rect(color = "black"))
gg_SES <- add_sensitivity_elements(gg_SES, cor_SES, legend_top=FALSE)
#loc_fig <- paste0(dir_figs,"BRS_SES_comparison")
#ggsave(paste0(loc_fig,".png"), gg_SES, width=100, height=100, units="mm", dpi=300)


## corr(SES, BRS) ####
cor_SES_BRS <- cor.test(pheno_SES$f738, pheno_SES$cov_bvr)
gg_SES_BRS <- ggplot(pheno, aes(x=SES, y=`BRS-ALL-cov_bvr`)) +
  geom_boxplot() +
  labs(x = 'SES (Household income)',
       y = 'Behavioral Risk Score (BRS)') +
  theme(axis.text.x = element_text(size=4),
        axis.title = element_text(size=6))
#loc_fig <- paste0(dir_figs,"BRS_SES_boxplot")
#ggsave(paste0(loc_fig,".png"), width=180, height=120, units="mm", dpi=300)

### combined plot ####
ggarrange(plotlist = list(gg_SES_BRS, gg_SES), ncol=2, labels='AUTO',
          align='h', widths=c(1.25,1.75))
loc_fig <- paste0(dir_figs,"BRS_SES1")
ggsave(paste0(loc_fig,".png"), width=167, height=100, units="mm", dpi=300)
