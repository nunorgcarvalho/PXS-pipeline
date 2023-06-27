# Libraries and paths ####
library(tidyverse)
library(data.table)
source("paths.R")

# path to UKBB phenotype file with field 26285
loc_pheno3 <- "/n/groups/patel/uk_biobank/project_22881_671028/ukb671028.csv"

# PRS-PXS correlation ####

# extract PRS-T2D column
cols_keep <- c("eid","26285-0.0")
pheno3 <- as_tibble(fread(loc_pheno3, select=cols_keep)) %>%
  rename(PRS_T2D = `26285-0.0`)
# adds PRS-T2D to pheno
pheno <- as_tibble(fread(loc_pheno)) %>%
  left_join(pheno3, by=c("IID"="eid"))
pheno_T2D <- pheno %>% select("IID",contains("T2D")) %>% drop_na()

## Plot
cor1 <- cor.test(pheno_T2D$PRS_T2D, pheno_T2D$PXS_T2D)
text <- paste0("r==",round(cor1$estimate,4),"~~p<2.2%*%10^-16")
ggplot(pheno_T2D, aes(x=PRS_T2D,y=PXS_T2D)) +
  geom_point(shape=1, alpha=0.01) +
  geom_smooth(method="lm") +
  theme_light() +
  labs(x = "Polygenic Risk Score (PRS) for T2D",
       y = "PolyeXposure Score (PXS) for T2D onset") +
  annotate("text", label = text, parse=TRUE,
           x = min(pheno_T2D$PRS_T2D) + 0.1, hjust = 0,
           y = max(pheno_T2D$PXS_T2D) - 0.1, vjust = 1,)
loc_out <- paste0(dir_results, "figures/PRS_vs_PXS.png")
ggsave(loc_out, width = 3000, height = 2000, units = "px")
