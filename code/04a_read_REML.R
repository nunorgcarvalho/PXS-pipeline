# reads heritability and gencorr script outputs

# loads directories and packages ####
library(tidyverse)
library(data.table)
source('code/00_paths.R')

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

# genetic correlation BRS vs CRFs comparisons ####
col_BRS <- 'BRS-ALL-cov_bvr' # main BRS term


fields <- as_tibble(fread('scratch/fields_tbl.txt')) %>%
  add_row(term='CRS-ALL', use_type='CRF',
          traitname = 'Clinical Risk Score (CRS)')
col_CRFs <- fields$term[fields$use_type == 'CRF']

CI95_z <- qnorm(0.975)
gc1 <- gencorr_REML %>%
  filter(trait1 == col_BRS,
         trait2 %in% col_CRFs) %>%
  left_join(fields %>% filter(use_type=='CRF') %>%
              select(trait2=term, traitname), by='trait2') %>%
  arrange(-abs(rg)) %>%
  mutate(rg_low = rg - CI95_z*rg_err,
         rg_upp = rg + CI95_z*rg_err,
         sign = ifelse(rg>0,'Positive','Negative'),
         traitname = factor(traitname))

ggplot(gc1, aes(x=factor(traitname, levels=traitname), y=abs(rg))) +
  geom_col(aes(fill = sign)) +
  geom_errorbar(aes(ymin=abs(rg_low), ymax=abs(rg_upp)),
                width=0.5) +
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=20)) +
  scale_y_continuous(expand=c(0,0,0.1,0)) +
  labs(x = 'Clinical Risk Factors and Score',
       y = 'Absolute genetic correlation (95% CI) with Behavioral Risk Score (BRS)',
       fill = 'Sign') +
  theme_bw() +
  theme(legend.position = 'top')
  

# genetic correlation of CRS ####
gc2 <- gencorr_REML %>%
  filter(trait1 == 'CRS-ALL') %>%
  add_row(gencorr_REML %>% filter(trait2 == 'CRS-ALL') %>%
            mutate(trait2 = trait1, trait1 = 'CRS-ALL')) %>%
  mutate(term = ifelse(trait2 == col_BRS, trait2,
                       str_replace_all(trait2,'_','\\.'))) %>%
  left_join(fields %>% filter(use_type=='behavior') %>%
              mutate(traitname = str_remove(traitname, "\\s*\\([^)]*\\)")) %>%
              add_row(term = col_BRS, traitname = 'Behavioral Risk Score (BRS)') %>%
              select(term, traitname), by='term') %>%
  arrange(-abs(rg)) %>%
  mutate(rg_low = rg - CI95_z*rg_err,
         rg_upp = rg + CI95_z*rg_err,
         sign = ifelse(rg>0,'Positive','Negative'),
         traitname = factor(traitname))

ggplot(gc2, aes(x=factor(traitname, levels=traitname), y=abs(rg))) +
  geom_col(aes(fill = sign)) +
  geom_errorbar(aes(ymin=abs(rg_low), ymax=abs(rg_upp)),
                width=0.5) +
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=40)) +
  scale_y_continuous(expand=c(0,0,0.1,0)) +
  labs(x = 'Behaviors and BRS',
       y = 'Absolute genetic correlation (95% CI) with Clinical Risk Score (CRS)',
       fill = 'Sign') +
  theme_bw() +
  theme(legend.position = 'top',
        axis.text.y = element_text(size=5))


# Heritability of behaviors and score ####
col_bvrs <- fields$term[fields$use_type == 'behavior']
col_bvrs_ <- str_replace_all(col_bvrs,'\\.','_')
hg1 <- gencorr_REML %>%
  filter(trait2 %in% c(col_bvrs_,col_BRS)) %>%
  group_by(trait2) %>%
  summarize(h2 = mean(trait2_h2g),
            h2_err = mean(trait2_h2g_err)) %>%
  mutate(term = ifelse(trait2 == col_BRS, trait2,
                       str_replace_all(trait2,'_','\\.')),
         h2_low = h2 - CI95_z*h2_err,
         h2_upp = h2 + CI95_z*h2_err) %>%
  arrange(-h2) %>%
  left_join(fields %>% filter(use_type=='behavior') %>%
              mutate(traitname = str_remove(traitname, "\\s*\\([^)]*\\)")) %>%
              add_row(term = col_BRS, traitname = 'Behavioral Risk Score (BRS)') %>%
              select(term, traitname), by='term')


ggplot(hg1, aes(x=factor(traitname, levels=traitname), y=h2)) +
  geom_col(aes(fill = term==col_BRS)) +
  geom_errorbar(aes(ymin=abs(h2_low), ymax=abs(h2_upp)),
                width=0.5) +
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=40)) +
  scale_y_continuous(expand=c(0,0,0.1,0),
                     breaks = seq(0,1,5/100),
                     minor_breaks = seq(0,1,1/100)) +
  labs(x = 'Behaviors and BRS',
       y = 'Heritability (h2)') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.y = element_text(size=5))


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
