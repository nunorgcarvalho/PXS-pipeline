# Need at least 8GB of RAM to load full LMM data later

# reuses code from aggregate_results and geographical_clustering to create clean
# tables and figures for publication

#### Libraries and directories ####
library(tidyverse)
library(data.table)
library(gt)
library(webshot)
library(ggrepel)

dir_script <- "~/jobs/PXS_pipeline/code/"
#dir_scratch <- "~/scratch3/PXS_pipeline/"
dir_scratch <- "~/scratch3/08-01_PXS_pipeline/"
dir_data_showcase <- "~/scratch3/key_data/" # contains 'Data_Dictionary_Showcase.tsv' from UKBB

# general constants
rounding_decimals <- 3

#### GT theme ####
gt_theme <- function(data) {
  
  tab_options(
    data = data,
    table.font.name = "Helvetica"
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  )
}

#### Table 1 ####
REML_expo <- as_tibble(fread(paste0(dir_scratch,"REML_exposures_results.txt")))
LMM_expo <- as_tibble(fread(paste0(dir_scratch,"LMM_exposures_results.txt"))) %>%
  select(-fieldname)


table1 <- left_join(REML_expo, LMM_expo, by="field") %>%
  select(fieldname, h2g, h2g_err, lambda, N_above_bonferroni) %>%
  mutate(h2g_CI = paste0("(",round(h2g - 2*h2g_err,rounding_decimals),", ",
                         round(h2g + 2*h2g_err,rounding_decimals),")")) %>%
  arrange(-h2g) %>%
  drop_na() %>%
  mutate(n = row_number()) %>%
  select(n,fieldname, h2g, h2g_CI, lambda, N_above_bonferroni) %>%
  filter(row_number() <= 10)

table1gt <- gt(table1) %>%
  gt_theme() %>%
  # summary_rows(
  #   columns = c("h2g","lambda"),
  #   fns = list("Mean" = ~mean(.))
  # ) %>%
  cols_label(
    n = "#",
    fieldname = "Exposure",
    h2g = md("h<sup>2</sup>"),
    h2g_CI = md("h<sup>2</sup> 95% CI"),
    lambda = "Lambda",
    N_above_bonferroni = "# of Significant SNPs"
  ) %>%
  fmt_number(
    columns = c("h2g","lambda"),
    suffixing = F,
    decimals = rounding_decimals
  ) %>%
  cols_align(
    align = "left",
    columns = c("fieldname")
  ) %>%
  tab_header(
    title = "Heritability and Lambda Inflation Factor of Exposures",
    subtitle = paste0("Displaying top ",nrow(table1)," exposures with highest heritability")
  ) %>%
  tab_footnote(
    footnote = "P < 0.05 / 19400443 (number of total SNPs)",
    locations = cells_column_labels(columns = c("N_above_bonferroni"))
  )
  
table1gt
gtsave(table1gt, "table1.png","./figures/")

#### Table 2 ####
REML_PXS <- as_tibble(fread(paste0(dir_scratch,"REML_PXS_results.txt")))
LMM_PXS <- as_tibble(fread(paste0(dir_scratch,"LMM_PXS_results.txt"))) %>%
  select(-fieldname)


table2 <- left_join(REML_PXS, LMM_PXS, by="field") %>%
  select(fieldname, h2g, h2g_err, lambda, N_above_bonferroni) %>%
  mutate(h2g_CI = paste0("(",round(h2g - 2*h2g_err,rounding_decimals),", ",
                         round(h2g + 2*h2g_err,rounding_decimals),")")) %>%
  arrange(-h2g) %>%
  mutate(n = row_number()) %>%
  select(n,fieldname, h2g, h2g_CI, lambda, N_above_bonferroni) %>%
  drop_na() %>%
  filter(row_number() <= 12)

table2gt <- gt(table2) %>%
  gt_theme() %>%
  summary_rows(
    columns = c("h2g","lambda"),
    fns = list("Mean" = ~mean(.)),
    formatter = fmt_number,
    decimals = rounding_decimals
  ) %>%
  cols_label(
    n = "#",
    fieldname = "PXS",
    h2g = md("h<sup>2</sup>"),
    h2g_CI = md("h<sup>2</sup> 95% CI"),
    lambda = "Lambda",
    N_above_bonferroni = "# of Significant SNPs"
  ) %>%
  fmt_number(
    columns = c("h2g","lambda"),
    suffixing = F,
    decimals = rounding_decimals
  ) %>%
  cols_align(
    align = "left",
    columns = c("fieldname")
  ) %>%
  tab_header(
    title = "Heritability and Lambda Inflation Factor of Exposures",
    subtitle = "Displaying 12 cardiometabolic health outcome PXSs"
  ) %>%
  tab_footnote(
    footnote = "P < 0.05 / 19400443 (number of total SNPs)",
    locations = cells_column_labels(columns = c("N_above_bonferroni"))
  )

table2gt
gtsave(table2gt, "table2.png","./figures/")

#### Table 3 ####
genCorrs <- as_tibble(fread(paste0(dir_scratch,"genCorr_CRF_results.txt")))


table3 <- genCorrs %>%
  select(CRF_fieldname, h2g2, h2g2_err, gencorr, gencorr_err) %>%
  mutate(h2g2_CI = paste0("(",round(h2g2 - 2*h2g2_err,rounding_decimals),", ",
                         round(h2g2 + 2*h2g2_err,rounding_decimals),")"),
         gencorr_CI = paste0("(",round(gencorr - 2*gencorr_err,rounding_decimals),", ",
                          round(gencorr + 2*gencorr_err,rounding_decimals),")")) %>%
  arrange(-gencorr) %>%
  mutate(n = row_number()) %>%
  select(n,CRF_fieldname, h2g2, h2g2_CI, gencorr, gencorr_CI) %>%
  drop_na() %>%
  filter(row_number() <= 6)

table3gt <- gt(table3) %>%
  gt_theme() %>%
  summary_rows(
    columns = c("h2g2","gencorr"),
    fns = list("|Mean|" = ~mean(abs(.))),
    formatter = fmt_number,
    decimals = rounding_decimals
  ) %>%
  cols_label(
    n = "#",
    CRF_fieldname = "Clinical Risk Factor",
    h2g2 = md("h<sup>2</sup>"),
    h2g2_CI = md("h<sup>2</sup> 95% CI"),
    gencorr = md("r<sub>g</sub>"),
    gencorr_CI = md("r<sub>g</sub> 95% CI")
  ) %>%
  fmt_number(
    columns = c("h2g2","gencorr"),
    suffixing = F,
    decimals = rounding_decimals
  ) %>%
  cols_align(
    align = "left",
    columns = c("CRF_fieldname")
  ) %>%
  tab_header(
    title = "Genetic Correlation between PXS for T2D and Clinical Risk Factors",
    subtitle = "Displaying 6 Clinical Risk Factors of Type 2 Diabetes"
  )

table3gt
gtsave(table3gt, "table3.png","./figures/")

#### Figure 1 ####
