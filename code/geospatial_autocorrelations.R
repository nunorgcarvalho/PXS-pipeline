library(tidyverse)
library(data.table)

dir_script <- "~/jobs/PXS_pipeline/code/"
dir_scratch <- "~/scratch3/PXS_pipeline/"
#loc_pheno_full1 <- "/n/groups/patel/uk_biobank/main_data_34521/ukb34521.tab"
loc_pheno_full2 <- "/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv"

cols_to_keep2 <- c("eid",paste0(c(129,130,22702,22704),"-0.0"))
pheno2 <- as_tibble(fread(loc_pheno_full2, select = cols_to_keep2))

ggplot(pheno2, aes(x=`22702-0.0`,y=`22704-0.0`)) +
  geom_point(alpha=0.01)
