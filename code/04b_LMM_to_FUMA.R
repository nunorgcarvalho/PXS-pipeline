# filters, selects, and compresses BOLT-LMM output to make it easier to 
# submit to FUMA


# libraries and paths ####
library(tidyverse)
library(data.table)
source('code/00_paths.R')

dir_LMM <- paste0(dir_scratch, 'LMM/')
dir_out <- paste0(dir_LMM,'FUMA_format/')
dir.create(dir_out, showWarnings = FALSE)

# set here the regex pattern of files to look for in the LMM/ directory
#regex_pattern = "LMM\\.f.*\\.bgen\\.txt"
regex_pattern <- '.*f.*\\.bgen\\.txt'
LMM_files <- list.files(dir_LMM, pattern=regex_pattern)
LMM_files

# loops through each file ####
for (i in 1:length(LMM_files)) {
  file <- LMM_files[i]
  print(paste(i,file))
  raw <- as_tibble(fread(paste0(dir_LMM, file)))
  
  LMM_out <- raw %>%
    select(SNP, CHR, BP, A0=ALLELE0, A1=ALLELE1, BETA, SE, P=P_BOLT_LMM_INF) %>%
    filter(P <0.1)
  
  loc_out <- paste0(dir_out, file, '.fuma.gz')
  
  fwrite(LMM_out, file=loc_out, sep='\t', compress='gzip')
}
