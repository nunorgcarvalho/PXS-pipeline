library(tidyverse)
library(data.table)

# Loads the phenotype file created in generate_PXS_coeffs.R
loc.pheno <- 'scratch/mortality/pheno_EC.death.txt'
pheno <- as_tibble(fread(loc.pheno))
coeffs <- as_tibble(fread('scratch/mortality/PXS_coefficients.txt'))
col.covs <- (coeffs %>%
               filter(substring(term,1,1) != 'f') %>%
               arrange(term))$term

# loads the file containing PXS coefficients
coeffs.tbl <- as_tibble(fread('scratch/mortality/PXS_coefficients.txt'))


pheno.wide <- pheno[c("FID","IID",coeffs.tbl$term)] %>% drop_na()
# 371,834 remaining

coeffs.vec <- coeffs.tbl$estimate
num_cols_skip <- ncol(pheno.wide) - nrow(coeffs.tbl)
PXSs <- c()
# computes dot product of expanded phenotype coding and XWAS coefficients for
# each individual in pheno file
for (i in 1:nrow(pheno.wide)) {
  iid_pheno.vec <- pheno.wide[i,(num_cols_skip+1):ncol(pheno.wide)] %>%
    unlist(use.names=FALSE)
  
  PXS <- (coeffs.vec %*% iid_pheno.vec)[1]
  PXSs <- c(PXSs,PXS)
  
  if (i %% 5000 == 0) {
    print(paste0("Computed PXS for ",i/1000,"k individuals"))
  }
}
# saves shortened version of pheno table with just IIDs and PXS
col.PXS <- 'PXS_death'
out.PXS <- pheno.wide %>% select(FID,IID)
#out.PXS[col.PXS] <- PXSs
out.PXS[col.PXS] <- (PXSs - mean(PXSs)) / sd(PXSs)


loc.out <- 'scratch/mortality/PXS_death.table.txt'
fwrite(out.PXS,loc.out,sep=" ")

# appends PXS to phenoEC & fullT2D tables
pheno <- pheno %>% left_join(out.PXS, by=c("FID","IID"))
fwrite(pheno, loc.pheno, sep='\t')

# saves NA IID
IIDs_NAs <- pheno[is.na(pheno[[col.PXS]]),c("FID","IID")]
fwrite(IIDs_NAs, 'scratch/mortality/IIDs_NA_exposures.txt', sep=' ', col.names = FALSE)
