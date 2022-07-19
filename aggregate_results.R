## Libraries and directories ##
library(tidyverse)
library(data.table)

dir_scratch <- "~/scratch3/PXS_pipeline/"
dir_script <- "~/jobs/PXS_pipeline/"
dir_data_showcase <- "~/scratch3/key_data/"

## Functions ##
extract_from_REML <- function(REML) {
  N <- length(REML)
  
  # extract time for analysis
  line <- REML[N]
  elapsed_hours <- as.numeric(str_split(line," ")[[1]][7]) / (60 * 60)
  
  # extract h2e
  line <- REML[N-4]
  line_split <- str_split(line, ": ")[[1]][2]
  h2e <- as.numeric(str_split(line_split, " ")[[1]][1])
  h2e_err <- parse_number(str_split(line_split, " ")[[1]][2])
  
  # extract h2g
  line <- REML[N-2]
  line_split <- str_split(line, ": ")[[1]][2]
  h2g <- as.numeric(str_split(line_split, " ")[[1]][1])
  h2g_err <- parse_number(str_split(line_split, " ")[[1]][2])
  
  out <- c(h2e, h2e_err, h2g, h2g_err, elapsed_hours)
  out
}

## Code ##

### REML results for each phenotypes
loc_phenolist <- paste0(dir_script,"phenotypes_ALL.txt")
pheno_list <- readLines(loc_phenolist)

REML_PXS_tbl <- tibble(
  field = as.character(),
  h2e = as.numeric(),
  h2e_err = as.numeric(),
  h2g = as.numeric(),
  h2g_err = as.numeric(),
  elapsed_hours = as.numeric()
)
REML_expo_tbl <- REML_PXS_tbl
for (pheno in pheno_list) {
  loc_REML <- paste0(dir_scratch,pheno,"/",pheno,"_PXS_BOLTREML.out")
  REML <- readLines(loc_REML)
  out <- extract_from_REML(REML)
  
  # adds row to table
  REML_PXS_tbl <- REML_PXS_tbl %>%
    add_row(
      field = pheno,
      h2e = out[1],
      h2e_err = out[2],
      h2g = out[3],
      h2g_err = out[4],
      elapsed_hours = out[5]
    )
  
  print(paste("Read REML results for",pheno))
}

### REML results for each exposure
loc_expolist <- paste0(dir_script,"exposures.txt")
exposures_list <- readLines(loc_expolist)

for (expo in exposures_list) {
  loc_REML <- paste0(dir_scratch,"exposures/",expo,"/",expo,"_BOLTREML.out")
  if (!file.exists(loc_REML)) {next}
  REML <- readLines(loc_REML)
  out <- extract_from_REML(REML)
  
  # adds row to table
  REML_expo_tbl <- REML_expo_tbl %>%
    add_row(
      field = expo,
      h2e = out[1],
      h2e_err = out[2],
      h2g = out[3],
      h2g_err = out[4],
      elapsed_hours = out[5]
    )
  
  print(paste("Read REML results for",expo))
}

### Appends fields name to REML tables
ukb_dict <- as_tibble(fread(paste0(dir_data_showcase,"Data_Dictionary_Showcase.tsv"))) %>%
  mutate(field = paste0("f.",FieldID,".0.0")) %>%
  select(field, fieldname=Field) %>%
  add_row(
    field = c("AF", "CAD", "COPD", "T2D","fev1_inst1"),
    fieldname = c("Atrial fibrillation", "Coronary artery disease", "Chronic obstructive pulmonary disease",
                  "Type 2 diabetes", "Forced expiratory volume in 1-second (FEV1)")
  )
REML_PXS_tbl <- REML_PXS_tbl %>% left_join(ukb_dict, by="field")
REML_expo_tbl <- REML_expo_tbl %>% left_join(ukb_dict, by="field")
