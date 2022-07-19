## Libraries and directories ##
library(tidyverse)
library(data.table)

dir_out <- "~/scratch3/PXS_pipeline/"
loc_coeffs <- "/n/groups/patel/yixuan/PXS_multi/all/coeffs_V9.csv"
loc_tally <- "/n/groups/patel/yixuan/PXS_multi/all/tally_X_V9.csv"
loc_pheno_full <- "/n/groups/patel/uk_biobank/main_data_9512/ukb9512.tab"
loc_pheno2_full <- "/n/groups/patel/uk_biobank/main_data_32144/ukb32144.csv"
loc_40PCs <- "~/scratch3/key_data/UKB_40PCs_500k.txt"
loc_CFRs_tbl <- "~/jobs/PXS_pipeline/CRFs_table.txt"

# List of exposures to remove from data for any reason, as well as covariates
remove_exposures <- c("f.20118.0.0")
covars <- c("f.31.0.0", # sex
            "f.34.0.0", # year of birth
            "f.54.0.0", # assessment center
            paste0("pc",1:40) # 40 PCs not from UKB field directly
            )
## Code

# loads coefficients and tally file from XWAS results
coeffs <- as_tibble(fread(loc_coeffs))
tally <- as_tibble(fread(loc_tally))
# loads list of CRFs
CRFs_tbl <- as_tibble(fread(loc_CFRs_tbl))

categorical_vars <- coeffs %>% select(X) %>%
  filter(X %in% tally$X) %>%
  group_by(X) %>% summarize(n_codings=n()) %>%
  mutate(
    data_type = ifelse(n_codings>1,"categorical","quantitative"),
    X = as.numeric(X)
  )

# create empty fields tibble
fields <- tibble(
  code = tally$X,
  field = paste0("f.",code,".0.0"),
  use_type = "exposure"
)

# Adds the covariates, exposures, CFRs, and IID fields to the fields tibble
fields <- fields %>%
  add_row(code=NA,field=covars,use_type="covar") %>%
  add_row(code=NA,field=CRFs_tbl$field_raw,use_type="CRF") %>%
  left_join(tally, by=c("code"="X")) %>%
  add_row(code = NA, field=c("f.eid","FID","IID"), use_type="identification",category="", fieldname="") %>%
  left_join(categorical_vars, by=c("code"="X")) %>%
  arrange(desc(field=="f.eid")) %>%
  arrange(desc(field=="IID")) %>%
  arrange(desc(field=="FID")) %>%
  filter(!(field %in% remove_exposures))

PCs40 <- as_tibble(fread(loc_40PCs)) %>% select(-FID)
pheno <- as_tibble(fread(loc_pheno_full, select=fields$field))
remaining_fields <- fields$field[!((fields$field) %in% colnames(pheno))]
pheno2 <- as_tibble(fread(loc_pheno2_full, select=c("eid",remaining_fields)))
# corrects field names to match others in pheno table, fields table, and CRFs table
for (i in 2:length(colnames(pheno2))) {
  field_raw <- colnames(pheno2[i])
  field <- CRFs_tbl[CRFs_tbl$field_raw==field_raw,"field"][[1]]
  colnames(pheno2)[i] <- field
  
  fields[fields$field==field_raw,"field"] <- field
}
pheno <- pheno %>%
  select(FID=f.eid, IID=f.eid,all_of(colnames(pheno[2:length(colnames(pheno))]))) %>%
  arrange(IID) %>%
  left_join(PCs40, by="IID") %>%
  left_join(pheno2, by=c("IID"="eid"))
  
# calculates the percentage of NA results per field and then removes fields with
# 100% missingness rate
for (i in 1:nrow(fields)) {
  field <- fields$field[i]
  
  if (!(field %in% colnames(pheno))) {pheno_NApct <- 1}
  else {pheno_NApct <- sum(is.na(pheno[field])) / nrow(pheno)}
  fields$pheno_NApct[i] <- pheno_NApct
}
fields <- fields %>% filter(pheno_NApct < 1)

# Saves phenotype file containing exposures and covariates
loc_out <- paste0(dir_out,"pheno_EC.txt")
write.table(pheno,loc_out,sep=" ",row.names=FALSE,quote=FALSE)

# Saves table with field information
loc_out <- paste0(dir_out,"fields_tbl.txt")
write.table(fields,loc_out,sep="\t",row.names=FALSE,quote=FALSE)
