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
dir.create('scratch/general_results/', showWarnings = FALSE)
fwrite(gencorr_REML, 'scratch/general_results/gencorr_REML.tsv', sep='\t')