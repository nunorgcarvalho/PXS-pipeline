# makes plots whose data comes primarily from the summary table

# loads directories and packages ####
library(tidyverse)
library(data.table)
library(ggrepel)
library(ggpubr)
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

# rg w/ CRS ####
gc2 <- sumtbl %>%
  select(term, shortname, rg=rg_CRS, rg_se=rg_se_CRS, beta) %>%
  arrange(-abs(rg)) %>%
  mutate(rg_low = rg - CI95_z*rg_se,
         rg_upp = rg + CI95_z*rg_se,
         sign = ifelse(rg>0,'Positive','Negative'),
         shortname = factor(shortname),
         rg_low_abs = ifelse(rg_low < 0 & rg>0, -Inf,abs(rg_low)),
         rg_upp_abs = ifelse(rg_upp > 0 & rg<0, -Inf,abs(rg_upp)),
         concordant = sign(rg) == sign(beta))

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

# rg w/ CRS vs T2D ####
gc5 <- gc2 %>%
  select(term, shortname, rg_CRS=rg, rg_CRS_low=rg_low, rg_CRS_upp=rg_upp) %>%
  left_join(gc4 %>% select(term=term2, rg_T2D=rg, rg_T2D_low=rg_low, rg_T2D_upp=rg_upp) %>%
              mutate(term = ifelse(term == col_BRS, 'BRS',term)),
            by='term') %>%
  left_join(sumtbl %>% select(term, rg_CRS_ldsc), by='term')

gc5_cor <- cor.test(gc5$rg_T2D, gc5$rg_CRS)
gc5_cor_ldsc <- cor.test(gc5$rg_T2D, gc5$rg_CRS_ldsc)
gc5_cor_CRS <- cor.test(gc5$rg_CRS, gc5$rg_CRS_ldsc)
gc5_label <- paste0('r = ', round(gc5_cor$estimate,3),
                    ' (p = ', formatC(gc5_cor$p.value,3),')')

gg_h2g_rgs_theme <- theme(
  axis.title = element_text(size=7),
  axis.text = element_text(size=6)
)

gg_gc5 <- ggplot(gc5, aes(x=rg_T2D, y=rg_CRS)) +
  geom_point(size=1) +
  geom_errorbarh(aes(xmin=rg_T2D_low, xmax=rg_T2D_upp), alpha=0.25) +
  geom_errorbar(aes(ymin=rg_CRS_low, ymax=rg_CRS_upp), alpha=0.25) +
  geom_text_repel(aes(label = shortname), size=1.5, seed = 2025,
                  force = 1, force_pull=1) +
  geom_vline(xintercept = 0, linetype='dashed', color=dash_color) +
  geom_hline(yintercept = 0, linetype='dashed', color=dash_color) +
  annotate('text', label = annotate_cor.test(gc5_cor), # throws out harmless warning message
           x=Inf, y=-Inf, hjust=1.05, vjust=-0.5, parse=FALSE, size=3) +
  labs(x = 'Genetic correlation (95% CI) with Type 2 Diabetes',
       y = 'Genetic correlation (95% CI) with Clinical Risk Score (CRS)') +
  theme_bw() + gg_h2g_rgs_theme

# h2 of bvrs/BRS ####
hg1 <- sumtbl %>%
  select(term, shortname, starts_with('h2')) %>%
  mutate(h2_low = h2 - CI95_z*h2_se,
         h2_upp = h2 + CI95_z*h2_se,
         h2_low_ldsc = h2_ldsc - CI95_z*h2_se_ldsc,
         h2_upp_ldsc = h2_ldsc + CI95_z*h2_se_ldsc) %>%
  arrange(-h2)
boldings <- hg1$term == 'BRS'
gg_hg1 <- ggplot(hg1, aes(x=factor(shortname, levels=shortname), y=h2)) +
  geom_col() +
  geom_errorbar(aes(ymin=abs(h2_low), ymax=abs(h2_upp)),
                width=0.5) +
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width=40)) +
  scale_y_continuous(expand=c(0,0,0.1,0),
                     breaks = seq(0,1,5/100),
                     minor_breaks = seq(0,1,1/100)) +
  labs(x = 'Behaviors and Behavioral Risk Score',
       y = 'SNP Heritability (95% CI)') +
  theme_bw() + gg_h2g_rgs_theme +
  theme(legend.position = 'none',
        axis.text.y = element_text(size=5,
                                   face = ifelse(boldings, 'bold','plain')))

# combined plot ####
ggarrange(plotlist = list(gg_hg1, gg_gc5), ncol=2, labels='AUTO',
          align='h', widths=c(1.5,2))
loc_fig <- paste0(dir_figs,"bvrs_h2g_rgs")
ggsave(paste0(loc_fig,".png"), width=240, height=140, units="mm", dpi=300)

rg_tbl %>% filter(term1 == 'CRS', term2 %in% c('BRS','T2D')) %>%
  select(term1,term2, rg, rg_se) %>%
  mutate(rg_low = rg - CI95_z*rg_se,
         rg_upp = rg + CI95_z*rg_se)


# BOLT-REML vs ldsc h2 estimates ####
lm_h2 <- lm(h2_ldsc ~ h2, data=hg1)
cor_h2 <- cor.test(hg1$h2, hg1$h2_ldsc)
ggplot(hg1, aes(x = h2, y=h2_ldsc)) +
  geom_abline(slope=1, color=dash_color, linetype='dashed') +
  geom_point() +
  geom_text_repel(aes(label = shortname), size=2) +
  geom_errorbarh(aes(xmin = h2_low, xmax = h2_upp), alpha=0.5) +
  geom_errorbar(aes(ymin = h2_low_ldsc, ymax = h2_upp_ldsc), alpha=0.5) +
  annotate('text', label = annotate_cor.test(cor_h2), # throws out harmless warning message
           x=Inf, y=-Inf, vjust=-0.5, hjust=1.05, parse=FALSE, size=3) +
  labs(x = 'h2 estimated by BOLT-REML (95% CI)',
       y = 'h2 estimated by BOLT-LDsc (95% CI)')
loc_fig <- paste0(dir_figs,"h2g_bvrs")
ggsave(paste0(loc_fig,".png"), width=160, height=160, units="mm", dpi=300)


# performance vs genetic metrics ####
gc3 <- sumtbl %>%
  select(term, shortname, beta_norm, AUC, AUC_se, rg_CRS, rg_CRS_se = rg_se_CRS,
         rg_T2D, rg_T2D_se = rg_se_T2D, h2, h2_se) %>%
  mutate(beta_norm_se = 0) %>%
  filter(term !='BRS')

# variable terms and labels
var1s <- list('beta_norm' = 'Effect on the BRS (normalized)',
              'AUC' = 'AUC for T2D')
var2s <- list('rg_CRS' = 'Genetic Correlation with CRS',
              'rg_T2D' = 'Genetic Correlation with T2D',
              'h2' = 'SNP Heritability')
gc3_plots <- list()
for (i in 1:length(var1s)) {
  var1 <- names(var1s)[i]
  for (j in 1:length(var2s)) {
    var2 <- names(var2s)[j]
    gc3_ij <- gc3 %>% select(term, shortname)
    gc3_ij$var1 <- gc3[[var1]]
    gc3_ij$var2 <- gc3[[var2]]
    # determines if absolute value is used for a variable
    var1_abs <- FALSE
    var2_abs <- FALSE
    if (var1 == 'AUC') { var2_abs = TRUE }
    if (var2 == 'h2' ) { var1_abs = TRUE }
    
    # sets absolute values
    if (var1_abs) {gc3_ij$var1 <- abs( gc3_ij$var1 )}
    if (var2_abs) {gc3_ij$var2 <- abs( gc3_ij$var2 )}
    
    # sets 95% CIs
    gc3_ij$var1_low <- gc3_ij$var1 - CI95_z * gc3[[paste0(var1,'_se')]]
    gc3_ij$var1_upp <- gc3_ij$var1 + CI95_z * gc3[[paste0(var1,'_se')]]
    gc3_ij$var2_low <- gc3_ij$var2 - CI95_z * gc3[[paste0(var2,'_se')]]
    gc3_ij$var2_upp <- gc3_ij$var2 + CI95_z * gc3[[paste0(var2,'_se')]]
    
    # gets correlation object
    cor_obj <- cor.test(gc3_ij$var1, gc3_ij$var2)
    # main plot
    gg0 <- ggplot(gc3_ij, aes(x=var2, y=var1)) +
      geom_smooth(method='lm', formula='y~x', se=FALSE) +
      geom_point(size=1) +
      geom_errorbar(aes(ymin = var1_low, ymax = var1_upp), alpha=0.25) +
      geom_errorbarh(aes(xmin = var2_low, xmax = var2_upp), alpha=0.25) +
      geom_text_repel(aes(label=shortname), max.overlaps = 5, size=2) +
      annotate('text', label = annotate_cor.test(cor_obj), # throws out harmless warning message
               x=-Inf, y=Inf, vjust=1.5, hjust=-0.05, parse=FALSE, size=3) +
      labs(x=var2s[j], y=var1s[i]) +
      theme_bw() +
      theme(axis.title = element_text(size=8),
            axis.text = element_text(size=6))
      
      
    # adds bars at x=0 or y=0
    if (!(var1_abs | var1=='AUC')) {gg0 <- gg0 + geom_hline(yintercept = 0, color=dash_color, linetype='dashed')}
    if (!(var2_abs | var2=='h2' )) {gg0 <- gg0 + geom_vline(xintercept = 0, color=dash_color, linetype='dashed')}
    
    # adds top/right ticks
    if (i == 1) { gg0 <- gg0 + scale_x_continuous(position = 'top')}
    if (j == 3) { gg0 <- gg0 + scale_y_continuous(position = 'right')}
    
    # removes y-axis labels for middle column
    if (!(j %in% c(1,3))) {
      gg0 <- gg0 + theme(axis.title.y = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks.y = element_blank())
    }
    
    # adds the word "Absolute" to axis label
    if (var1 == 'beta_norm' & var2 == 'h2') {gg0 <- gg0 + labs(y=paste0('Absolute ', var1s[i]))}
    if (var1 == 'AUC' & var2 != 'h2') {gg0 <- gg0 + labs(x=paste0('Absolute ', var2s[j]))}
    
    gc3_plots[[(i-1)*3 + j]] <- gg0
  }
}
ggarrange(plotlist=gc3_plots,ncol=3, nrow=2, labels = 'AUTO',
          hjust= 0.5 * c(-1,1,1, -1,1,1),
          vjust= c(rep(1.25,3), rep(0.25,3)),
          widths=c(1,0.9,1), heights=c(1,1))
loc_fig <- paste0(dir_figs,"beta_AUC_rg_h2")
ggsave(paste0(loc_fig,".png"), width=270, height=150, units="mm", dpi=300)


# AUC vs Effect Size ####
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