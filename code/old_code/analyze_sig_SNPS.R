## Libraries and directories ##
library(tidyverse)
library(data.table)

dir_script <- "~/jobs/PXS_pipeline/code/"
dir_scratch <- "~/scratch3/PXS_pipeline/"
dir_data_showcase <- "~/scratch3/key_data/"

## Code ##

loc_phenolist <- paste0(dir_script,"../input_data/phenotypes.txt")
pheno_list <- readLines(loc_phenolist)

# builds dictionary of exposures per phenotype
expos <- vector(mode="list", length=length(pheno_list))
names(expos) <- pheno_list
for (pheno in pheno_list) {
  loc_expos <- paste0(dir_scratch,pheno,"/",pheno,"_exposures.txt")
  exposures <- readLines(loc_expos)
  expos[[pheno]] <- exposures
}

bonferroni <- 0.05/ 19400443
sig_SNPs_PXS <- as_tibble(fread(paste0(dir_scratch, "LMM_PXS_sig_SNPs.txt"))) #%>% filter(P < bonferroni)
sig_SNPs_expo <- as_tibble(fread(paste0(dir_scratch, "LMM_expo_sig_SNPs.txt"))) #%>% filter(P < bonferroni)


PXS_expo_SNPs_agreement <- tibble(
  field = character(),
  n_snps_original = numeric(),
  n_snps_LD = numeric(),
  n_snps_found = numeric()
)
bin_size <- 100000
for (pheno in pheno_list) {
  print(paste("Looking at", pheno))
  exposures <- expos[[pheno]]
  sig_SNPs_pheno <- sig_SNPs_PXS %>%
    filter(field == pheno) %>%
    arrange(P) %>%
    mutate(BP_lower = BP - bin_size/2,
           BP_upper = BP + bin_size/2 - 1)
  
  n_snps_original <- nrow(sig_SNPs_pheno)
  
  sig_SNPs_LD <- sig_SNPs_pheno[0,]
  check <- TRUE
  while (nrow(sig_SNPs_pheno) > 0) {
    sig_SNPs_LD <- sig_SNPs_LD %>% add_row(sig_SNPs_pheno[1,])
    skip <- nrow(sig_SNPs_LD)
    lower <- sig_SNPs_LD$BP_lower[[skip]]
    upper <- sig_SNPs_LD$BP_upper[[skip]]
    chrom <- sig_SNPs_LD$CHR[[skip]]
    
    sig_SNPs_pheno <- sig_SNPs_pheno %>%
      filter( !( (BP >= lower) & (BP <= upper) & (CHR==chrom) ) )
  }
  n_snps_LD <- nrow(sig_SNPs_LD)
  
  sig_SNPs_expo_filter <- sig_SNPs_expo %>%
    filter(field %in% exposures)
  
  sig_SNPs_LD$contains_expo_SNP <- FALSE
  for (i in 1:nrow(sig_SNPs_expo_filter)) {
    js <- 1:nrow(sig_SNPs_LD)
    BP <- sig_SNPs_expo_filter$BP[[i]]
    CHR <- sig_SNPs_expo_filter$CHR[[i]]
    within_LD <- js[ which( (BP >= sig_SNPs_LD$BP_lower) & (BP <= sig_SNPs_LD$BP_upper) & (CHR == sig_SNPs_LD$CHR) ) ]
    sig_SNPs_LD$contains_expo_SNP[within_LD] <- TRUE
    sig_SNPs_expo_filter$within_LD[i] <- (length(within_LD) > 0)
    if (i %% 5000 == 0) {print(paste("Done with",i, "rows out of", nrow(sig_SNPs_expo_filter)))}
  }
  n_snps_found <- nrow(sig_SNPs_LD %>% filter(contains_expo_SNP))
  
  PXS_expo_SNPs_agreement <- PXS_expo_SNPs_agreement %>% add_row(
    field = pheno,
    n_snps_original = n_snps_original,
    n_snps_LD = n_snps_LD,
    n_snps_found = n_snps_found
  )
  
  print(paste0(n_snps_found," out of ",n_snps_LD," independent significant SNPs in PXS_",pheno," were found significant in the exposures"))
}
PXS_expo_SNPs_agreement <- PXS_expo_SNPs_agreement %>%
  mutate(portion_found = n_snps_found / n_snps_LD) %>%
  arrange(-portion_found)
