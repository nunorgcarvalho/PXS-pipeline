## Libraries and directories ##
library(tibble)
library(magrittr)
library(dplyr)
library(data.table)

loc_phenolist <- "phenotypes_ALL.txt"
loc_CRFs_tbl <- "CRFs_table.txt"
dir_scratch <- "~/scratch3/PXS_pipeline/"
dir_out <- "./"

## Code ###

# Part 1 - make_exposures_list
loc_fields <- paste0(dir_scratch,"fields_tbl.txt")

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

# Part 2 - prepare_PXS_CRF
CRFs_tbl <- as_tibble(fread(loc_CRFs_tbl))
loc_pheno <- paste0(dir_scratch,"pheno_EC.txt")
pheno_tbl <- as_tibble(fread(loc_pheno))
for (pheno in phenolist) {
  print(paste0("Appending PXS_",pheno," to pheno table"))
  
  # appends PXS
  loc_PXS <- paste0(dir_scratch,pheno,"/","PXS_",pheno,".txt")
  PXSs <- as_tibble(fread(loc_PXS)) %>% select(-FID)
  pheno_tbl <- pheno_tbl %>% left_join(PXSs, by="IID")
  
  dir_pheno <- paste0(dir_scratch,pheno,"/")
  # writes list of exposures
  exposures <- (fields %>%
                  filter(use_type=="exposure") %>%
                  select(field, all_of(pheno)) %>%
                  mutate(sum = rowSums(across(all_of(pheno)))) %>%
                  filter(sum > 0) )$field
  
  loc_out <- paste0(dir_pheno,pheno,"_exposures.txt")
  file_out <- file(loc_out)
  writeLines(exposures,file_out)
  close(file_out)
  
  # writes list of CRFs
  CRFs_pheno <- CRFs_tbl %>% filter(phenotype==pheno)
  if (nrow(CRFs_pheno)==0) {next}
  CRFs <- CRFs_pheno$field
  loc_out <- paste0(dir_pheno,pheno,"_CRFs.txt")
  file_out <- file(loc_out)
  writeLines(CRFs,file_out)
  close(file_out)
}
write.table(pheno_tbl, loc_pheno, sep=" ", quote=FALSE, row.names=FALSE)
