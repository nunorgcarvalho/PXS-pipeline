# reads heritability and gencorr script outputs

# loads directories and packages ####
library(tidyverse)
library(data.table)
source('code/00_paths.R')

# general files and variables ####
shortnames <- as_tibble(fread('input_data/term_shortnames.tsv'))
fields <- as_tibble(fread('scratch/fields_tbl.txt'))
col_CRFs <- fields$term[fields$use_type == 'CRF']
col_bvrs <- fields$term[fields$use_type == 'behavior']
col_bvrs_ <- str_replace_all(col_bvrs,'\\.','_')
col_BRS <- 'BRS-ALL-cov_bvr' # main BRS term
CI95_z <- qnorm(1 - (1 - 0.95)/2)

# REML gencorrs ####
dir_gencorrs <- paste0(dir_scratch, 'gencorrs/')
rg_files_out <- list.files(dir_gencorrs, pattern='^REML.*\\.out$')


gencorr_REML <- tibble(file='', h2_cohort='',
                       trait1='',trait2='',
                       trait1_h2g=0, trait1_h2g_err=0,
                       trait2_h2g=0, trait2_h2g_err=0,
                       rg=0, rg_err=0, elapsed_hours=0)[0,]
for (rg_file in rg_files_out) {
  print(rg_file)
  file_split <- str_split(rg_file,'\\.')[[1]]
  if (length(file_split) == 5) {
    h2_cohort <- file_split[2]
    file_split <- file_split[-2]
  } else {h2_cohort <- NA}
  rg_out <- readLines(paste0(dir_gencorrs, rg_file))
  
  ## reads rg out file ####
  N <- length(rg_out)
  
  # extract time for analysis
  line <- rg_out[N]
  elapsed_hours <- as.numeric(str_split(line," ")[[1]][7]) / (60 * 60)
  
  # extract h2g1
  line <- rg_out[N-4]
  line_split <- str_split(line, ": ")[[1]][2]
  h2g1 <- as.numeric(str_split(line_split, " ")[[1]][1])
  h2g1_err <- parse_number(str_split(line_split, " ")[[1]][2])
  
  # extract rg_out
  line <- rg_out[N-3]
  line_split <- str_split(line, ": ")[[1]][2]
  gencorr <- as.numeric(str_split(line_split, " ")[[1]][1])
  gencorr_err <- parse_number(str_split(line_split, " ")[[1]][2])
  
  # extract h2g2
  line <- rg_out[N-2]
  line_split <- str_split(line, ": ")[[1]][2]
  h2g2 <- as.numeric(str_split(line_split, " ")[[1]][1])
  h2g2_err <- parse_number(str_split(line_split, " ")[[1]][2])
  
  
  gencorr_REML <- gencorr_REML %>% add_row(
    file=rg_file, h2_cohort=h2_cohort,
    trait1=file_split[2], trait2=file_split[3],
    trait1_h2g=h2g1, trait1_h2g_err=h2g1_err,
    trait2_h2g=h2g2, trait2_h2g_err=h2g2_err,
    rg=gencorr, rg_err=gencorr_err, elapsed_hours=elapsed_hours
  )
}

# behavior summary table ####

# extracts significant behaviors for BRS-ALL-cov_bvr

## coefficients of terms for BRS ####
load('scratch/BRS_models/cv_glm_models.RData')
BRS_cvglm <- models[[str_remove(col_BRS, 'BRS-')]]
BRS_coeffs <- as_tibble(fread('scratch/BRS_models/BRS_coefficients.txt')) %>%
  select(term, traitname, beta = all_of(paste0('beta-',col_BRS))) %>%
  drop_na() %>% 
  mutate(factor_SD = BRS_cvglm$factor_SDs,
         beta_norm = beta * factor_SD,
         beta_norm_abs = abs(beta_norm),
         effect_sign = sign(beta_norm)) %>%
  filter(beta != 0)


bvr_sumtbl <- tibble(term = col_bvrs[col_bvrs %in% BRS_coeffs$term]) %>%
  left_join(shortnames, by='term') %>%
  left_join(fields %>% select(term, data_type), by='term') %>%
  left_join(BRS_coeffs %>% select(term, beta, beta_norm), by='term') %>%
  left_join(gencorr_REML %>% 
              mutate(term = str_replace_all(trait2,'_','\\.')) %>%
              filter(term %in% BRS_coeffs$term, trait1 == 'CRS-ALL') %>%
              select(term, h2=trait2_h2g, h2_err = trait2_h2g_err,
                     rg_CRS = rg, rg_err_CRS = rg_err))

dir.create('scratch/general_results/', showWarnings = FALSE)
fwrite(bvr_sumtbl, 'scratch/general_results/bvr_sumtbl.txt', sep='\t')


# genetic correlation BRS vs CRFs comparisons ####
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
  theme(legend.position = 'top')
  

# genetic correlation with CRS ####
gc2 <- bvr_sumtbl %>%
  select(term, shortname, rg=rg_CRS, rg_err=rg_err_CRS, beta) %>%
  add_row(gencorr_REML %>% filter(trait2=='CRS-ALL') %>%
            select(term=trait1, rg, rg_err) %>%
            mutate(beta=1) %>%
            left_join(shortnames %>% select(-traitname), by='term')) %>%
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
                                   face = ifelse(boldings, 'bold','plain')))


# Heritability of behaviors and score ####

hg1 <- bvr_sumtbl %>%
  select(term, shortname, h2, h2_err) %>%
  add_row(gencorr_REML %>% filter(trait2==col_BRS) %>%
            rename(term=trait2) %>% group_by(term) %>%
            summarize(h2 = mean(trait2_h2g),
                      h2_err = mean(trait2_h2g_err)) %>%
            left_join(shortnames %>% select(-traitname), by='term')) %>%
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

# coefficients ####
library(ggrepel)

gc3 <- bvr_sumtbl %>%
  rename(rg=rg_CRS, rg_err=rg_err_CRS) %>%
  mutate(rg_low = rg - CI95_z*rg_err,
         rg_upp = rg + CI95_z*rg_err)

gc3_cor <- cor.test(gc3$rg, gc3$beta_norm)
ggplot(gc3, aes(x=beta_norm, y=rg)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  #geom_smooth(method='lm', se=FALSE) +
  geom_point() +
  geom_errorbar(aes(ymin = rg_low, ymax = rg_upp), alpha=0.5) +
  geom_text_repel(aes(label = shortname), size=2) +
  theme_bw() +
  labs(x = 'Effect on BRS (normalized)',
       y = 'Genetic correlation (95%) w/ Clinical Risk Score (CRS)',
       subtitle = paste0('r = ', round(gc3_cor$estimate,3),
                         ' :: p = ', formatC(gc3_cor$p.value,3)))

#
hg2 <- bvr_sumtbl %>%
  mutate(h2_low = h2 - CI95_z*h2_err,
         h2_upp = h2 + CI95_z*h2_err)

hg2_cor <- cor.test(hg2$h2, abs(hg2$beta_norm))
ggplot(hg2, aes(x=abs(beta_norm), y=h2)) +
  geom_point() +
  geom_errorbar(aes(ymin = h2_low, ymax = h2_upp), alpha=0.5) +
  geom_text_repel(aes(label = shortname), size=2) +
  theme_bw() +
  labs(x = 'Absolute effect on BRS (normalized)',
       y = 'Heritability (95%) w/',
       subtitle = paste0('r = ', round(hg2_cor$estimate,3),
                         ' :: p = ', formatC(hg2_cor$p.value,3)))

hg2b_cor <- cor.test(hg2$h2, log10(abs(hg2$beta_norm)) )
ggplot(hg2, aes(x=log10(abs(beta_norm)), y=h2)) +
  geom_point() +
  geom_errorbar(aes(ymin = h2_low, ymax = h2_upp), alpha=0.5) +
  geom_text_repel(aes(label = shortname), size=2) +
  theme_bw() +
  labs(x = 'Log10 of Absolute effect on BRS (normalized)',
       y = 'Heritability (95%) w/',
       subtitle = paste0('r = ', round(hg2b_cor$estimate,3),
                         ' :: p = ', formatC(hg2b_cor$p.value,3)))



#############################################################
# clean table ####
Cstats <- as_tibble(fread('scratch/BRS_models/BRS_Cstats.txt'))

## h2 of different BRSs####
h2_table <- bind_rows(
  gencorr_REML %>% select(h2_cohort, trait=trait1, h2g=trait1_h2g, h2g_err=trait1_h2g_err),
  gencorr_REML %>% select(h2_cohort, trait=trait2, h2g=trait2_h2g, h2g_err=trait2_h2g_err)
) %>% #drop_na() %>%
  mutate(training_cohort = sapply(str_split(trait, "-"), `[`, 2),
         model_factors = sapply(str_split(trait, "-"), `[`, 3)) %>%
  left_join(Cstats %>% select(training_cohort, model_factors, model_label) %>% distinct(),
            by=c('training_cohort','model_factors')) %>%
  select(h2_cohort, training_cohort, model_label, h2g, h2g_err) %>%
  group_by(h2_cohort, training_cohort, model_label) %>%
  summarize(h2g = mean(h2g), h2g_err = mean(h2g_err))

ggplot(h2_table, aes(x=model_label, color=training_cohort)) +
  geom_point(aes(y=h2g), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = h2g - 2*h2g_err, ymax = h2g + 2*h2g_err),
                position = position_dodge(width = 0.5), width=0.5) +
  facet_wrap(~ h2_cohort, nrow=1) +
  theme_bw() +
  theme(legend.position = 'top') +
  labs(x = 'Model Factors',
       y = 'h2 (95% CI)',
       color = 'BRS trained in:',
       subtitle = 'Panels denote the cohort that h2 was estimated from')



## gencorr ####
rg_table <- gencorr_REML %>%
  mutate(training_cohort = sapply(str_split(trait1, "-"), `[`, 2),
         trait1_model_factors = sapply(str_split(trait1, "-"), `[`, 3),
         trait2_model_factors = sapply(str_split(trait2, "-"), `[`, 3)) %>%
  left_join(Cstats %>%
              select(training_cohort, trait1_model_factors = model_factors,
                     trait1_model_label = model_label) %>% distinct(),
            by=c('training_cohort','trait1_model_factors')) %>%
  left_join(Cstats %>%
              select(training_cohort, trait2_model_factors = model_factors,
                     trait2_model_label = model_label) %>% distinct(),
            by=c('training_cohort','trait2_model_factors')) %>%
  select(h2_cohort, training_cohort, trait1_model_label, trait2_model_label, rg, rg_err)
