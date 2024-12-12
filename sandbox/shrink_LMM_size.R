library(tidyverse)
library(data.table)

# shrinks the size of LMM files so when Gzipped, they're <600Mb
file.path <- 'scratch/T2D/LMM_PXS_T2D_BMIadj_bgen.txt'
LMM.raw <- as_tibble(fread(file.path))

LMM.shrunk <- LMM.raw %>%
  select(-c(GENPOS, INFO, CHISQ_LINREG,P_LINREG))

file.out.path <- paste0(file.path,'.shrunk')
fwrite(LMM.shrunk, file.out.path, sep='\t')
