# Convert a large BOLT-LMM summary stats file to something easily input by FUMA (gzipped)
# FUMA: https://fuma.ctglab.nl/

# library(getopt)
library(tidyverse)
library(data.table)
PVAL_THRESHOLD <- 0.5


# spec = matrix(c(
#   'file_in', 'i', 1, "character",
#   'sex_file_in', 's', 1, "character",
#   'path_out', 'p', 1, "character"
# ), byrow=TRUE, ncol=4)
# opt = getopt(spec)

#file_in <- opt$file_in
file_in <- "/n/groups/patel/nuno/PXS-pipeline/scratch/T2D/LMM_T2D_all_bgen.txt"
#sex_file_in <- opt$sex_file_in
sex_file_in <- NULL
#path_out <- opt$path_out
path_out <- "/n/groups/patel/nuno/PXS-pipeline/scratch/T2D/"

print(file_in)
print(path_out)
print(basename(file_in))
#file_out <- sprintf('%s.fuma.gz', basename(file_in))
file_out <- paste0(path_out, basename(file_in), ".fuma.gz")
print(file_out)

gwas_stats <- as_tibble(fread(file_in)) %>%
  select(SNP, CHR, BP, A0=ALLELE0, A1=ALLELE1, P=P_BOLT_LMM_INF, Beta=BETA, SE) %>%
  filter(P < PVAL_THRESHOLD)
#gwas_stats_to_fuma <- gwas_stats %>% select(SNP, CHR, BP, ALLELE0, ALLELE1, P_BOLT_LMM_INF, BETA, SE) %>% rename(A1 = ALLELE1, A0 = ALLELE0, P=P_BOLT_LMM_INF, Beta=BETA) %>% filter(P < PVAL_THRESHOLD)

if(!is.null(sex_file_in)) {
  gwas_stats <- rbind(gwas_stats, read_delim(sex_file_in, delim = "\t") %>% select(SNP, CHR, BP, ALLELE0, ALLELE1, P_BOLT_LMM_INF, BETA, SE) %>% rename(A1 = ALLELE1, A0 = ALLELE0, P=P_BOLT_LMM_INF, Beta=BETA) %>% filter(P < PVAL_THRESHOLD))
}

fwrite(gwas_stats, file=file_out, sep="\t", compress = "gzip", nThread = 10, scipen=50)
