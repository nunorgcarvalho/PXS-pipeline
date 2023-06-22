# Libraries and paths ####
library(tidyverse)
library(data.table)
library(gt)

dir_results <- paste0(dir_script,"../final_results/")

# Common theme(s)
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
rounding_decimals <- 3

# GWAS Catalog Table ####
## Trait count ####
loc_tbl <- paste0(dir_results, "GWAS_catalog_trait_count.txt")
tbl_raw <- as_tibble(fread(loc_tbl))
tbl <- tbl_raw %>% arrange(-n) %>% filter(row_number() <= 25)

tbl_gt <- gt(tbl) %>%
  gt_theme() %>%
  cols_label(Trait = "GWAS Catalog Trait", n = "# of SNP replications") %>%
  cols_align(align = "left", columns = c("Trait")) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = everything())
  )
tbl_gt
gtsave(tbl_gt, "GWAS_catalog_trait_count.png",dir_results)

## Group count ####
loc_tbl <- paste0(dir_results, "GWAS_catalog_group_count.txt")
tbl_raw <- as_tibble(fread(loc_tbl))
tbl <- tbl_raw %>% arrange(-n) %>% filter(row_number() <= 25)

tbl_gt <- gt(tbl) %>%
  gt_theme() %>%
  cols_label(Group = "GWAS Catalog Trait Group", n = "# of SNP replications") %>%
  cols_align(align = "left", columns = c("Group")) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = everything())
  )
tbl_gt
gtsave(tbl_gt, "GWAS_catalog_group_count.png",dir_results)
