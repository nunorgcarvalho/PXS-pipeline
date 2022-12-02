## Libraries and directories ##
library(tibble)
library(magrittr)
library(dplyr)
library(data.table)

dir_script <- "~/jobs/PXS_pipeline/code/"
dir_scratch <- "~/scratch3/PXS_pipeline/"

## Code ###

# Part 1 - make_exposures_list
loc_fields <- paste0(dir_scratch,"fields_tbl.txt")
loc_coeffs <- paste0(dir_script,"../input_data/PXS_coefficients.txt")

loc_phenolist <- paste0(dir_script,"../input_data/phenotypes.txt")
phenolist <- readLines(loc_phenolist)
fields <- as_tibble(fread(loc_fields))
coeffs <- as_tibble(fread(loc_coeffs))
# exposures <- (fields %>%
#   filter(use_type=="exposure") %>%
#   select(field, all_of(phenolist)) %>%
#   mutate(sum = rowSums(across(all_of(phenolist)))) %>%
#   filter(sum > 0) )$field
exposures <- (coeffs %>%
  select(term, disease) %>%
  left_join(fields %>% select(term,use_type), by="term") %>%
  filter(use_type=="exposure", disease %in% phenolist))$term

loc_out <- paste0(dir_script,"../input_data/exposures.txt")
file_out <- file(loc_out)
writeLines(exposures,file_out)
close(file_out)

# Part 2 - prepare_PXS_CRF
loc_CRFs_tbl <- paste0(dir_script,"../input_data/CRFs_table.txt")
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
  exposures <- (coeffs %>%
                  select(term, disease) %>%
                  left_join(fields %>% select(term,use_type), by="term") %>%
                  filter(use_type=="exposure", disease == pheno))$term
  
  loc_out <- paste0(dir_pheno,pheno,"_exposures.txt")
  file_out <- file(loc_out)
  writeLines(exposures,file_out)
  close(file_out)
  
  # makes list of IIDs with missing PXS
  IIDs_NAs <- pheno_tbl %>% select(FID, IID, PXS = all_of(paste0("PXS_",pheno))) %>%
                 filter(is.na(PXS)) %>% select(FID,IID)
  loc_out <- paste0(dir_pheno,"IIDs_NA_exposures.txt")
  write.table(IIDs_NAs, loc_out, row.names = FALSE, col.names = FALSE, quote = FALSE)
  #file_out <- file(loc_out)
  #writeLines(IIDs_NAs,file_out)
  #close(file_out)
  
  # writes list of CRFs
  CRFs_pheno <- CRFs_tbl %>% filter(phenotype==pheno)
  if (nrow(CRFs_pheno)==0) {next}
  CRFs <- paste0("f",sapply(str_split(CRFs_pheno$field,"\\."), '[', 2))
  loc_out <- paste0(dir_pheno,pheno,"_CRFs.txt")
  file_out <- file(loc_out)
  writeLines(CRFs,file_out)
  close(file_out)
}
write.table(pheno_tbl, loc_pheno, sep=" ", quote=FALSE, row.names=FALSE)
