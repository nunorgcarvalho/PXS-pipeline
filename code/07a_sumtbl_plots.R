# makes plots whose data comes primarily from the summary table

# loads directories and packages ####
library(tidyverse)
library(data.table)
library(ggrepel)
source('code/00_paths.R')
source('code/00_plotting.R')


# main data sources ####
sumtbl <- as_tibble(fread('scratch/general_results/sumtbl.tsv'))
fields <- as_tibble(fread('scratch/fields_tbl.txt'))
gencorr_REML <- as_tibble(fread('scratch/general_results/gencorr_REML.tsv'))
rg_tbl <- as_tibble(fread(paste0(dir_scratch,'general_results/behavior_rg_table.tsv')))
shortnames <- as_tibble(fread('input_data/term_shortnames.tsv'))

# shared variables ####
col_CRFs <- fields$term[fields$use_type == 'CRF']
col_BRS <- 'BRS-ALL-cov_bvr' # main BRS term
CI95_z <- qnorm(1 - (1 - 0.95)/2)

# rg BRS vs CRFs comparisons ####
gc1 <- gencorr_REML %>%
  filter(trait1 == col_BRS,
         trait2 %in% col_CRFs) %>%
  left_join(shortnames %>% select(trait2=term, shortname), by='trait2') %>%
  arrange(-abs(rg)) %>%
  mutate(rg_low = rg - CI95_z*rg_err,
         rg_upp = rg + CI95_z*rg_err,
         sign = ifelse(rg>0,'Positive','Negative'),
         shortname = factor(shortname))

ggplot(gc1, aes(x=factor(shortname, levels=shortname), y=abs(rg))) +
  geom_col(aes(fill = sign)) +
  geom_errorbar(aes(ymin=abs(rg_low), ymax=abs(rg_upp)),
                width=0.5) +
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width=20)) +
  scale_y_continuous(expand=c(0,0,0.1,0)) +
  labs(x = 'Clinical Risk Factors and Score',
       y = 'Absolute genetic correlation (95% CI) with Behavioral Risk Score (BRS)',
       fill = 'Sign') +
  theme_bw() +
  theme(legend.position = 'top',
        legend.margin = margin(b = -5))
loc_fig <- paste0(dir_figs,"gencorr_BRS_CRFs")
ggsave(paste0(loc_fig,".png"), width=180, height=120, units="mm", dpi=300)


# rg w/ CRS ####
gc2 <- sumtbl %>%
  select(term, shortname, rg=rg_CRS, rg_err=rg_err_CRS, beta) %>%
  arrange(-abs(rg)) %>%
  mutate(rg_low = rg - CI95_z*rg_err,
         rg_upp = rg + CI95_z*rg_err,
         sign = ifelse(rg>0,'Positive','Negative'),
         shortname = factor(shortname),
         rg_low_abs = ifelse(rg_low < 0 & rg>0, -Inf,abs(rg_low)),
         rg_upp_abs = ifelse(rg_upp > 0 & rg<0, -Inf,abs(rg_upp)),
         concordant = sign(rg) == sign(beta))

# vector denoting which traitnames to bold
boldings <- gc2$term == col_BRS

ggplot(gc2, aes(x=factor(shortname, levels=shortname), y=abs(rg))) +
  geom_col(aes(fill = concordant)) +
  geom_errorbar(aes(ymin=rg_low_abs, ymax=rg_upp_abs),
                width=0.5) +
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width=40)) +
  scale_y_continuous(expand=c(0,0,0.1,0)) +
  labs(x = 'Behaviors and BRS',
       y = 'Absolute genetic correlation (95% CI) with Clinical Risk Score (CRS)',
       fill = 'Genetic correlation concordant with BRS effect') +
  theme_bw() +
  theme(legend.position = 'top',
        axis.text.y = element_text(size=6,
                                   face = ifelse(boldings, 'bold','plain')),
        legend.margin = margin(b = -5))

loc_fig <- paste0(dir_figs,"gencorr_CRS_bvrs")
ggsave(paste0(loc_fig,".png"), width=180, height=120, units="mm", dpi=300)


# rg w/ T2D ####
gc4 <- rg_tbl %>% filter(term1 == 'T2D' | term2 == 'T2D',
                         term1 != 'CRS') %>%
  select(-starts_with('h2g'),-starts_with('intercept'), -file) %>%
  mutate(term2 = str_replace_all(term2,'_','\\.')) %>%
  mutate(term2 = ifelse(term2=='BRS',col_BRS,term2)) %>%
  left_join(shortnames, by=c('term2'='term')) %>%
  left_join(sumtbl %>% select(term,beta_norm), by=c('term2'='term')) %>%
  arrange(-abs(rg)) %>%
  mutate(rg_low = rg - CI95_z*rg_se,
         rg_upp = rg + CI95_z*rg_se,
         sign = ifelse(rg>0,'Positive','Negative'),
         shortname = factor(shortname),
         rg_low_abs = ifelse(rg_low < 0 & rg>0, -Inf,abs(rg_low)),
         rg_upp_abs = ifelse(rg_upp > 0 & rg<0, -Inf,abs(rg_upp)),
         concordant = sign(rg) == sign(beta_norm))

ggplot(gc4, aes(x=factor(shortname, levels=shortname), y=abs(rg))) +
  geom_col(aes(fill = concordant)) +
  geom_errorbar(aes(ymin=rg_low_abs, ymax=rg_upp_abs), width=0.5) +
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width=40)) +
  scale_y_continuous(expand=c(0,0,0.1,0)) +
  labs(x = 'Behaviors and BRS',
       y = 'Absolute genetic correlation (95% CI) with Type 2 Diabetes',
       fill = 'Genetic correlation concordant with BRS effect') +
  theme_bw() +
  theme(legend.position = 'top',
        axis.text.y = element_text(size=6,
                                   face = ifelse(boldings, 'bold','plain')),
        legend.margin = margin(b = -5))
loc_fig <- paste0(dir_figs,"gencorr_T2D_bvrs")
ggsave(paste0(loc_fig,".png"), width=180, height=120, units="mm", dpi=300)

# rg w/ CRS vs T2D ####
gc5 <- gc2 %>%
  select(term, shortname, rg_CRS=rg, rg_CRS_low=rg_low, rg_CRS_upp=rg_upp) %>%
  left_join(gc4 %>% select(term=term2, rg_T2D=rg, rg_T2D_low=rg_low, rg_T2D_upp=rg_upp),
            by='term')

gc5_cor <- cor.test(gc5$rg_T2D, gc5$rg_CRS)
gc5_label <- paste0('r = ', round(gc5_cor$estimate,3),
                    ' (p = ', formatC(gc5_cor$p.value,3),')')
ggplot(gc5, aes(x=rg_T2D, y=rg_CRS)) +
  geom_point() +
  geom_errorbarh(aes(xmin=rg_T2D_low, xmax=rg_T2D_upp), alpha=0.25) +
  geom_errorbar(aes(ymin=rg_CRS_low, ymax=rg_CRS_upp), alpha=0.25) +
  geom_text_repel(aes(label = shortname), size=2) +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_hline(yintercept = 0, linetype='dashed') +
  annotate('text', label = gc5_label,
           x = Inf, y=-Inf, hjust=1.05, vjust=-1) +
  labs(x = 'Genetic correlation (95% CI) with Type 2 Diabetes',
       y = 'Genetic correlation (95% CI) with Clinical Risk Score (CRS)') +
  theme_bw()
#loc_fig <- paste0(dir_figs,"gencorr_T2D_vs_CRS")
#ggsave(paste0(loc_fig,".png"), width=180, height=180, units="mm", dpi=300)


# h2 of bvrs/BRS ####
hg1 <- sumtbl %>%
  select(term, shortname, h2, h2_err) %>%
  mutate(h2_low = h2 - CI95_z*h2_err,
         h2_upp = h2 + CI95_z*h2_err) %>%
  arrange(-h2)

ggplot(hg1, aes(x=factor(shortname, levels=shortname), y=h2)) +
  geom_col() +
  geom_errorbar(aes(ymin=abs(h2_low), ymax=abs(h2_upp)),
                width=0.5) +
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width=40)) +
  scale_y_continuous(expand=c(0,0,0.1,0),
                     breaks = seq(0,1,5/100),
                     minor_breaks = seq(0,1,1/100)) +
  labs(x = 'Behaviors and BRS',
       y = 'Heritability (h2)') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.y = element_text(size=6,
                                   face = ifelse(boldings, 'bold','plain')))
loc_fig <- paste0(dir_figs,"h2g_bvrs")
ggsave(paste0(loc_fig,".png"), width=180, height=180, units="mm", dpi=300)

# effect vs rg w/ CRS ####
gc3 <- sumtbl %>%
  filter(term != col_BRS) %>%
  rename(rg=rg_CRS, rg_err=rg_err_CRS) %>%
  mutate(rg_low = rg - CI95_z*rg_err,
         rg_upp = rg + CI95_z*rg_err)

gc3_cor <- cor.test(gc3$rg, gc3$beta_norm)
ggplot(gc3, aes(x=beta_norm, y=rg)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_point() +
  geom_errorbar(aes(ymin = rg_low, ymax = rg_upp), alpha=0.5) +
  geom_text_repel(aes(label = shortname), size=2) +
  theme_bw() +
  labs(x = 'Effect on BRS (normalized)',
       y = 'Genetic correlation (95%) w/ Clinical Risk Score (CRS)',
       subtitle = paste0('r = ', round(gc3_cor$estimate,3),
                         ' :: p = ', formatC(gc3_cor$p.value,3)))
loc_fig <- paste0(dir_figs,"gencorr_CRS_vs_effect")
ggsave(paste0(loc_fig,".png"), width=180, height=180, units="mm", dpi=300)

# effect vs h2 ####
hg2 <- sumtbl %>%
  filter(term != col_BRS) %>%
  mutate(h2_low = h2 - CI95_z*h2_err,
         h2_upp = h2 + CI95_z*h2_err)

hg2_cor <- cor.test(hg2$h2, abs(hg2$beta_norm))
ggplot(hg2, aes(x=abs(beta_norm), y=h2)) +
  geom_point() +
  geom_errorbar(aes(ymin = h2_low, ymax = h2_upp), alpha=0.5) +
  geom_text_repel(aes(label = shortname), size=2) +
  theme_bw() +
  labs(x = 'Absolute effect on BRS (normalized)',
       y = 'Heritability (95%)',
       subtitle = paste0('r = ', round(hg2_cor$estimate,3),
                         ' :: p = ', formatC(hg2_cor$p.value,3)))

# log-transformed effect
hg2b_cor <- cor.test(hg2$h2, log10(abs(hg2$beta_norm)) )
ggplot(hg2, aes(x=log10(abs(beta_norm)), y=h2)) +
  geom_point() +
  geom_errorbar(aes(ymin = h2_low, ymax = h2_upp), alpha=0.5) +
  geom_text_repel(aes(label = shortname), size=2) +
  theme_bw() +
  labs(x = 'Log10 of Absolute effect on BRS (normalized)',
       y = 'Heritability (95%)',
       subtitle = paste0('r = ', round(hg2b_cor$estimate,3),
                         ' :: p = ', formatC(hg2b_cor$p.value,3)))


# AUC vs Effect Size
auc1 <- sumtbl %>%
  filter(term != 'BRS')
ROC_tbl <- as_tibble(fread(paste0(dir_scratch,'general_results/ROC_tbl.tsv'))) %>%
  filter(term %in% c('cov','BRS'))

auc1_cor <- cor.test(auc1$AUC, abs(auc1$beta_norm))
ggplot(auc1, aes(x=abs(beta_norm), y=AUC)) +
  annotate('rect', fill='black', alpha=0.1, xmin=-Inf, xmax=Inf,
           ymin = ROC_tbl$AUC_low[1], ymax = ROC_tbl$AUC_upp[1]) +
  geom_hline(yintercept=ROC_tbl$AUC[1], linetype='dashed') + # cov-only
  annotate('rect', fill='black', alpha=0.1, xmin=-Inf, xmax=Inf,
           ymin = ROC_tbl$AUC_low[2], ymax = ROC_tbl$AUC_upp[2]) +
  geom_hline(yintercept=ROC_tbl$AUC[2], linetype='dashed') + # BRS
  geom_point() +
  geom_errorbar(aes(ymin = AUC_low, ymax = AUC_upp), alpha=0.5) +
  geom_text_repel(aes(label = shortname), size=2) +
  annotate('text', x=Inf, y=ROC_tbl$AUC[1], vjust=-0.5, hjust=1.05, size=3,
           label='AUC for covariate-only model') +
  annotate('text', x=Inf, y=ROC_tbl$AUC[2], vjust=1.5, hjust=1.05, size=3,
           label='AUC for BRS model') +
  annotate('text', label = annotate_cor.test(auc1_cor), # throws out harmless warning message
           x=-Inf, y=ROC_tbl$AUC_low[2], vjust=1.05, hjust=-0.05, parse=FALSE) +
  theme_bw() +
  labs(x = 'Absolute effect on BRS (normalized)',
       y = 'AUC for Type 2 Diabetes')
loc_fig <- paste0(dir_figs,"bvrs_AUC_vs_effect")
ggsave(paste0(loc_fig,".png"), width=180, height=150, units="mm", dpi=300)