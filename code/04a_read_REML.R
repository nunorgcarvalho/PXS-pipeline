# reads heritability and gencorr script outputs

# loads directories and packages ####
library(tidyverse)
library(data.table)
source('code/00_paths.R')

# REML gencorrs ####
dir_gencorrs <- paste0(dir_scratch, 'gencorrs/')
rg_files_out <- list.files(dir_gencorrs, pattern='*.out')

gencorr_REML <- tibble(file='', h2_cohort='',
                       trait1='',trait2='',
                       trait1_h2g=0, trait1_h2g_err=0,
                       trait2_h2g=0, trait2_h2g_err=0,
                       rg=0, rg_err=0, elapsed_hours=0)[0,]
for (rg_file in rg_files_out) {
  print(rg_file)
  file_split <- str_split(rg_file,'\\.')[[1]]
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
    file=rg_file, h2_cohort=file_split[2],
    trait1=file_split[3], trait2=file_split[4],
    trait1_h2g=h2g1, trait1_h2g_err=h2g1_err,
    trait2_h2g=h2g2, trait2_h2g_err=h2g2_err,
    rg=gencorr, rg_err=gencorr_err, elapsed_hours=elapsed_hours
  )
}


# clean table ####
Cstats <- as_tibble(fread('scratch/BRS_models/BRS_Cstats.txt'))

## h2 ####
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
