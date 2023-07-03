# Libraries and paths ####
library(tidyverse)
library(data.table)
library(gt)
source('code/paths.R')

dir_poster <- paste0(dir_results,"poster/")

# Common code ####
gt_theme <- function(data) {
  tab_options(
    data = data,
    table.font.name = "Helvetica"
  ) %>%
    cols_align(
      align = "center",
      columns = everything()
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels(columns = everything())
    ) %>%
    sub_missing(missing_text = "")
}
rnd_dec <- 3
shortnames <- as_tibble(fread(paste0(dir_script,"../input_data/exposures_shortnames.csv")))

# Heritabilities + Lambdas + # of genetic loci ####
h2_tbl <- as_tibble(fread(paste0(dir_results,"h2_ldsc_REML.txt")))
Genes_tbl <- as_tibble(fread(paste0(dir_results, "gene_associations_n.txt")))
GWAS_tbl <- as_tibble(fread(paste0(dir_results,"GWAS_summary_tbl.txt")))

b1 <- qnorm(1 - (0.025 / nrow(h2_tbl)))
tbl <- h2_tbl %>%
  left_join(Genes_tbl, by="shortname") %>%
  left_join(GWAS_tbl, by="term") %>%
  select(shortname, REML_h2, n_sig_SNPs, n) %>%
  arrange(-REML_h2)#arrange(shortname)

tbl_gt <- gt(tbl) %>%
  gt_theme() %>%
  cols_label(
    shortname = "Behavior",
    REML_h2 = md("h<sup>2</sup>"),
    n_sig_SNPs = "Significant SNPs",
    n = "Significant Genes",
  ) %>%
  fmt_number(columns = c("REML_h2"), decimals = 3) %>%
  cols_width(n ~ px(80),
             n_sig_SNPs ~ px(80)) %>%
  cols_align(align = "left", columns = c("shortname"))
tbl_gt
gtsave(tbl_gt, "behaviors_genetic_profile.png",dir_poster)
## Statistics
(tbl %>% filter(shortname != "PXS for T2D onset"))$REML_h2 %>% mean()

# Genetic correlations ####
loc_expolist <- paste0(dir_script,"../input_data/exposures.txt")
exposures_list <- readLines(loc_expolist)

REML_results <- as_tibble(fread(paste0(dir_scratch, "genCorr_REML_results.txt")))

REML_results %>% filter(pheno1_term %in% exposures_list) %>%
  group_by(pheno2_term,pheno2_shortname) %>%
  summarize(h2 = mean(h2g2, na.rm=TRUE),
            mean_abs_gencorr= mean(abs(gencorr), na.rm=TRUE)) %>%
  arrange(-mean_abs_gencorr)

REML_results %>% #filter(pheno1_term %in% exposures_list) %>%
  group_by(pheno1_term,pheno1_shortname) %>%
  summarize(h2 = mean(h2g1, na.rm=TRUE),
            mean_abs_gencorr= mean(abs(gencorr), na.rm=TRUE)) %>%
  arrange(-mean_abs_gencorr) %>% print(n=Inf)

# GTEx expression ####
GTEx_general <- as_tibble(fread(paste0(dir_results,"GTEx_general.txt")))
GTEx_general %>% group_by(VARIABLE) %>%
  filter(term != "PXS_T2D", significant_P_adj) %>%
  summarize(n = n()) %>%
  arrange(-n)

GTEx_specific <- as_tibble(fread(paste0(dir_results,"GTEx_specific.txt")))
GTEx_specific %>% group_by(FULL_NAME) %>%
  filter(term != "PXS_T2D", significant_P_adj) %>%
  summarize(n = n()) %>%
  arrange(-n)
