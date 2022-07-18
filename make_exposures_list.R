## Libraries and directories ##
library(tibble)
library(magrittr)
library(dplyr)
library(data.table)

loc_phenolist <- "phenotypes.txt"
loc_fields <- "~/scratch3/PXS_pipeline//fields_tbl.txt"
dir_out <- "./"

## Code ###

phenolist <- readLines(loc_phenolist)
fields <- tibble::as_tibble(fread(loc_fields))
exposures <- (fields %>%
  filter(use_type=="exposure") %>%
  select(field, all_of(phenolist)) %>%
  mutate(sum = rowSums(across(all_of(phenolist)))) %>%
  filter(sum > 0) )$field

loc_out <- paste0(dir_out,"exposures.txt")
file_out <- file(loc_out)
writeLines(exposures,file_out)
close(file_out)
