# Libraries and paths ####
library(tidyverse)
library(data.table)
library(gt)
source('code/paths.R')

dir_tbls <- paste0(dir_results,"tables/")

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

# GWAS Catalog Table ####
## Trait count ####
loc_tbl <- paste0(dir_results, "GWAS_catalog_trait_count.txt")
tbl_raw <- as_tibble(fread(loc_tbl))
tbl <- tbl_raw %>% arrange(-n) %>% filter(row_number() <= 25)
tbl_csv <- tbl %>% rename(`# of SNP replications`=n)
fwrite(tbl_csv, paste0(dir_tbls,"GWAS_catalog_trait_count.csv"))

# tbl_gt <- gt(tbl) %>%
#   gt_theme() %>%
#   cols_label(Trait = "GWAS Catalog Trait", n = "# of SNP replications") %>%
#   cols_align(align = "left", columns = c("Trait"))
# tbl_gt
# gtsave(tbl_gt, "GWAS_catalog_trait_count.png",dir_tbls)
# gtsave(tbl_gt, "GWAS_catalog_trait_count.rtf",dir_tbls)

## Group count ####
loc_tbl <- paste0(dir_results, "GWAS_catalog_group_count.txt")
tbl_raw <- as_tibble(fread(loc_tbl))
tbl <- tbl_raw %>% arrange(-n) %>% filter(row_number() <= 25)
tbl_csv <- tbl %>% rename(`# of SNP replications`=n)
fwrite(tbl_csv, paste0(dir_tbls,"GWAS_catalog_group_count.csv"))

# tbl_gt <- gt(tbl) %>%
#   gt_theme() %>%
#   cols_label(Group = "GWAS Catalog Trait Group", n = "# of SNP replications") %>%
#   cols_align(align = "left", columns = c("Group"))
# tbl_gt
# gtsave(tbl_gt, "GWAS_catalog_group_count.png",dir_tbls)
# gtsave(tbl_gt, "GWAS_catalog_group_count.rtf",dir_tbls)

# Fields Table ####
## All behavioral traits ####
loc_tbl <- paste0(dir_script,"../input_data/UKBB_behaviors.txt")
tbl_raw <- as_tibble(fread(loc_tbl))
ukb_dict <- as_tibble(fread(paste0(dir_data_showcase,"Data_Dictionary_Showcase.tsv")))

tbl <- tbl_raw %>% filter(Exposure_Class=="agency") %>%
  arrange(FieldID) %>%
  left_join(ukb_dict[,c("FieldID","ValueType")], by="FieldID") %>%
  select(FieldID, Field, ValueType)
tbl_csv <- tbl %>% rename(`UKBB Field ID` = FieldID,
                          `Field Name` = Field,
                          `Data Type` = ValueType)
fwrite(tbl_csv, paste0(dir_tbls,"fields_ALL_behaviors.csv"))

# tbl_gt <- gt(tbl) %>%
#   gt_theme() %>%
#   cols_label(
#     FieldID = "UKBB Field ID",
#     Field = "Field Name",
#     ValueType = "Data Type") %>%
#   cols_align(align = "left", columns = c("Field"))
# tbl_gt
# gtsave(tbl_gt, "fields_ALL_behaviors.png",dir_tbls) # maybe just save as .csv?
# gtsave(tbl_gt, "fields_ALL_behaviors.rtf",dir_tbls)

## PXS Coefficients ####
loc_tbl <- paste0(dir_script,"../input_data/PXS_coefficients.txt")
tbl_raw <- as_tibble(fread(loc_tbl))
fields <- as_tibble(fread(paste0(dir_scratch,"fields_tbl.txt")))

# sets maps for renaming covariates and use types
map_shortname <- setNames(c("Sex","Age",paste0("PC",1:40)),
                          c("sex","age",paste0("pc",1:40)) )
map_usetype <- setNames(c("Covariate", "Behavior"),
                        c("covar","exposure") )

# joins coefficients, shortnames, and field usetype
tbl <- tbl_raw %>%
  left_join(shortnames, by="term") %>%
  left_join(fields[c("term","use_type")], by="term") %>%
  mutate(shortname = ifelse(is.na(shortname),term,shortname)) %>%
  mutate(shortname = recode(shortname, !!!map_shortname),
         use_type = recode(use_type, !!!map_usetype)) %>%
  select(use_type, fieldID, value, shortname, estimate, p.value)
# manually sets field ID for sex and age
tbl$fieldID[tbl$shortname=="Sex"] <- 31
tbl$fieldID[tbl$shortname=="Age"] <- 34
# resorts table to first be exposures (sorted by p-value) and then covariates
tbl <-  tbl %>% filter(use_type=="Behavior") %>% arrange(p.value) %>%
  add_row(tbl %>%
            filter(use_type=="Covariate") %>%
            mutate(shortname = factor(shortname, levels = unname(map_shortname))) %>%
            arrange(shortname)
          )
tbl_csv <- tbl
fwrite(tbl_csv, paste0(dir_tbls,"PXS_coefficients.csv"))


# tbl_gt <- gt(tbl) %>%
#   gt_theme() %>%
#   cols_label(
#     use_type = "Variable Type",
#     fieldID = "UKBB Field ID",
#     value = "UKBB Value",
#     shortname = "Variable Name",
#     estimate = "Standardized Coefficient",
#     p.value = "XWAS p-value"
#   ) %>%
#   fmt_number(columns = c("estimate"), decimals = 5) %>%
#   fmt_scientific(columns = "p.value") %>%
#   cols_width(fieldID ~ px(80),
#              value ~ px(60),
#              estimate ~ px(100))
# tbl_gt
# gtsave(tbl_gt, "PXS_coefficients.png",dir_tbls) # maybe just save as .csv?
# gtsave(tbl_gt, "PXS_coefficients.rtf",dir_tbls)

# Heritabilities + Lambdas + # of genetic loci ####
h2_tbl <- as_tibble(fread(paste0(dir_results,"h2_ldsc_REML.txt")))
Genes_tbl <- as_tibble(fread(paste0(dir_results, "gene_associations_n.txt")))
Loci_tbl <- as_tibble(fread(paste0(dir_results, "genomic_risk_loci.txt"))) %>%
  group_by(shortname) %>% summarize(n_sig_loci = n()) %>% arrange(-n_sig_loci)
GWAS_tbl <- as_tibble(fread(paste0(dir_results,"GWAS_summary_tbl.txt")))

b1 <- qnorm(1 - (0.025 / nrow(h2_tbl)))
tbl <- h2_tbl %>%
  left_join(Genes_tbl, by="shortname") %>%
  left_join(Loci_tbl, by="shortname") %>%
  left_join(GWAS_tbl, by="term") %>%
  mutate(REML_h2_low = REML_h2 - b1 * REML_h2_err,
         REML_h2_upp = REML_h2 + b1 * REML_h2_err) %>%
  mutate(REML_h2_CI = paste0("[",sprintf(paste0("%.",rnd_dec,"f"), REML_h2_low),", ",
                             sprintf(paste0("%.",rnd_dec,"f"), REML_h2_upp), "]"),
         shortname = factor(shortname, levels=shortnames$shortname)) %>%
  select(shortname, REML_h2, REML_h2_CI, lambda_BOLT, n_sig_SNPs, n, n_sig_loci) %>%
  arrange(-REML_h2)#arrange(shortname)

tbl_gt <- gt(tbl) %>%
  gt_theme() %>%
  cols_label(
    shortname = "Behavior",
    REML_h2 = md("h<sup>2</sup>"),
    REML_h2_CI = md("h<sup>2</sup> 95% CI"),
    lambda_BOLT = "Lambda",
    n_sig_SNPs = "Significant SNPs",
    n = "Significant Genes",
    n_sig_loci = "Significant Loci"
  ) %>%
  fmt_number(columns = c("REML_h2","lambda_BOLT"), decimals = 3) %>%
  cols_width(n ~ px(120)) %>%
  cols_align(align = "left", columns = c("shortname"))
tbl_gt
gtsave(tbl_gt, "behaviors_genetic_profile.png",dir_tbls)
gtsave(tbl_gt, "behaviors_genetic_profile.rtf",dir_tbls)