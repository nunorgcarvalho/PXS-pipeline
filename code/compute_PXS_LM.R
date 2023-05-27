## Libraries and directories ##
disease <- commandArgs(trailingOnly = TRUE)[1]
source('paths.R')
dir_out <- paste0(dir_scratch,disease,"/")
setwd(dir_out)
sink(paste0(disease,"_compute_PXS_LM.Rout"))

library(tidyverse)
library(data.table)

## Code ###

#col_coeff <- paste0("estimate_",disease)
col_coeff <- "estimate"

# Loads the phenotype file created in generate_PXS_coeffs.R
pheno <- as_tibble(fread(loc_phenoEC))
loc_fields <- "../fields_tbl.txt"
fields <- as_tibble(fread(loc_fields)) 
col_covs <- (fields %>% filter(use_type=="covar"))$term

# loads the file containing XWAS coefficients and filters to just to the disease
loc_coeffs <- paste0(dir_script,"../input_data/PXS_coefficients.txt")
coeffs <- as_tibble(fread(loc_coeffs)) %>%
  select(term, estimate, disease) #%>%

var_tbl <- coeffs %>% left_join(fields, by="term")
# saves list of all exposures
exposures <- (var_tbl %>% filter(use_type=="exposure"))$term
loc_out <- paste0(dir_script,"../input_data/exposures.txt")
file_out <- file(loc_out)
writeLines(exposures,file_out)
close(file_out)

var_tbl <-var_tbl %>%
  filter(disease == .GlobalEnv$disease) %>%
  arrange(use_type, term)
# saves list of exposures specific to disease
exposures <- (var_tbl %>% filter(use_type=="exposure"))$term
loc_out <- paste0(dir_scratch,disease,"/",disease,"_exposures.txt")
file_out <- file(loc_out)
writeLines(exposures,file_out)
close(file_out)



pheno_wide <- pheno[c("FID","IID",var_tbl$term)] %>%
  drop_na() # removes NA values

coeffs_vec <- var_tbl$estimate
coeffs_vec[is.na(coeffs_vec)] <- 0
num_cols_skip <- ncol(pheno_wide) - nrow(var_tbl)
PXSs <- c()
# computes dot product of expanded phenotype coding and XWAS coefficients for
# each individual in pheno file
for (i in 1:nrow(pheno_wide)) {
  iid_pheno_vec <- pheno_wide[i,(num_cols_skip+1):ncol(pheno_wide)] %>%
    unlist(use.names=FALSE)
  
  PXS <- (coeffs_vec %*% iid_pheno_vec)[1]
  PXSs <- c(PXSs,PXS)
  
  if (i %% 5000 == 0) {
    print(paste0("Computed PXS for ",i/1000,"k individuals"))
  }
}
# saves shortened version of pheno table with just IIDs and PXS
col_PXS <- paste0("PXS_",disease)
out_PXS <- pheno_wide %>% select(FID,IID)
out_PXS[col_PXS] <- PXSs

# checks normality
# ggplot(out_PXS, aes(x=PXS_T2D)) +
#   geom_density( color="red") +
#   stat_function(fun = dnorm, n = 101, args = list(mean = mean(out_PXS$PXS_T2D), sd = sd(out_PXS$PXS_T2D)))
# forces standardized PXS value:
out_PXS[col_PXS] <- (PXSs - mean(PXSs)) / sd(PXSs)


loc_out <- paste0(dir_out,"PXS_",disease,".txt")
fwrite(out_PXS,loc_out,sep=" ")
# appends PXS to phenoEC & fullT2D tables
pheno <- pheno %>% left_join(out_PXS, by=c("FID","IID"))
fwrite(pheno, loc_phenoEC, sep="\t", na="NA", quote=FALSE)

loc_fullT2D <- paste0(dir_scratch, "phenoEC_fullT2D.txt")
T2D_definitions_out <- as_tibble(fread(loc_fullT2D))
T2D_definitions_out <- T2D_definitions_out %>% left_join(out_PXS, by=c("FID","IID"))
fwrite(T2D_definitions_out, loc_fullT2D, sep="\t", na="NA", quote=FALSE)

print("Done computing PXS")

# makes list of IIDs with missing PXS
if (col_PXS %in% colnames(pheno)) {
  pheno <- pheno %>% select(-all_of(col_PXS))
}
pheno <- pheno %>% left_join(out_PXS, by=c("FID","IID"))
IIDs_NAs <- pheno[is.na(pheno[[col_PXS]]),c("FID","IID")]
loc_out <- paste0(dir_scratch,disease,"/IIDs_NA_exposures.txt")
write.table(IIDs_NAs, loc_out, row.names = FALSE, col.names = FALSE, quote = FALSE)

### Computes env LM: PXS ~ sex + age + assessment_center + PCs
PXS_lm_tbl <- pheno %>% select(FID,IID,all_of(col_covs)) %>%
  right_join(out_PXS, by=c("FID","IID")) %>%
  rename("PXS" = all_of(col_PXS)) %>%
  mutate(sex = as.factor(sex),
         assessment_center = as.factor(assessment_center)) %>%
  select(-FID,-IID) %>%
  drop_na()

PXS_lm <- lm(data=PXS_lm_tbl, PXS ~ .)
loc_out <- paste0("PXS_",disease,"_envLM.rds")
saveRDS(PXS_lm,loc_out)
summary(PXS_lm)

print("Done computing envLM")

# writes list of CRFs
loc_CRFs_tbl <- paste0(dir_script,"../input_data/CRFs_table.txt")
CRFs_tbl <- as_tibble(fread(loc_CRFs_tbl))
CRFs_pheno <- CRFs_tbl %>% filter(phenotype==disease)
if (nrow(CRFs_pheno)==0) {next}
CRFs <- paste0("f",sapply(str_split(CRFs_pheno$field,"\\."), '[', 2))
loc_out <- paste0(dir_scratch,disease,"/",disease,"_CRFs.txt")
file_out <- file(loc_out)
writeLines(CRFs,file_out)
close(file_out)

