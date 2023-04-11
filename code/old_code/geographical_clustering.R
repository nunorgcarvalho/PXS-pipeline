## Libraries and directories ##
#library(plyr)
library(tidyverse)
library(data.table)
library(GGally)
# Functions #
source("~/jobs/PXS_pipeline/code/helper_functions.R")

dir_script <- "~/jobs/PXS_pipeline/code/"
#dir_scratch <- "~/scratch3/PXS_pipeline/"
dir_scratch <- "~/scratch3/08-01_PXS_pipeline/"
dir_data_showcase <- "~/scratch3/key_data/" # contains 'Data_Dictionary_Showcase.tsv' and 'Codings.tsv' from UKBB






## Code ##
ukb_dict <- as_tibble(fread(paste0(dir_data_showcase,"Data_Dictionary_Showcase.tsv"))) %>%
  mutate(field = paste0("f.",FieldID,".0.0")) %>%
  select(field, fieldname=Field) %>%
  add_row(
    field = c("AF", "CAD", "COPD", "T2D","fev1_inst1"),
    fieldname = c("Atrial fibrillation", "Coronary artery disease", "Chronic obstructive pulmonary disease",
                  "Type 2 diabetes", "Forced expiratory volume in 1-second (FEV1)")
  )

AC_tbl <- as_tibble(fread(paste0(dir_data_showcase,"Codings.tsv"),quote="")) %>%
  filter(Coding == 10)

loc_pheno <- paste0(dir_scratch,"pheno_EC.txt")
pheno <- as_tibble(fread(loc_pheno)) %>%
  dplyr::rename(sex = f.31.0.0,
         year_of_birth = f.34.0.0,
         assessment_center = f.54.0.0) %>%
  mutate(sex = as.factor(sex),
         assessment_center = factor(assessment_center,
                                    levels = AC_tbl$Value,
                                    labels = AC_tbl$Meaning))

loc_phenolist <- paste0(dir_script,"../input_data/phenotypes.txt")
pheno_list <- readLines(loc_phenolist)

# Loops through XWAS diseases
ANOVA_PXS_tbl <- tibble(
  field = as.character(),
  f_stat = as.numeric(),
  p_value = as.numeric()
)
ANOVA_expo_tbl <- ANOVA_PXS_tbl
AC_PXS_means <- tibble(assessment_center = levels(pheno$assessment_center))
AC_expo_means <- AC_PXS_means
for (field in pheno_list) {
  
  col_PXS <- paste0("PXS_",field)
  loc_out <- paste0(dir_scratch,field,"/","ANOVAbp_",field,".png")
  out_list <- calculate_ANOVA(col_PXS, loc_out, ANOVA_PXS_tbl, AC_PXS_means)
  ANOVA_PXS_tbl <- out_list[[1]]
  AC_PXS_means <- out_list[[2]]
  print(paste("Calculated ANOVA and saved plot for", field))
  
}

# Loops through exposures
loc_fields <- paste0(dir_scratch,"fields_tbl.txt")
fields <- as_tibble(fread(loc_fields))

# only performs analysis on quantitative exposures found in the data
exposures_list <- (fields %>% filter(use_type == "exposure",
                                     data_type == "quantitative"))$field
exposures_list <- exposures_list[exposures_list %in% colnames(pheno)]

for (expo in exposures_list) {
  dir_expo <-  paste0(dir_scratch,"exposures/",expo)
  dir.create(dir_expo, recursive = TRUE)
  loc_out <- paste0(dir_expo,"/","ANOVAbp_",expo,".png")
  out_list <- calculate_ANOVA(expo, loc_out, ANOVA_expo_tbl, AC_expo_means, TRUE)
  ANOVA_expo_tbl <- out_list[[1]]
  AC_expo_means <- out_list[[2]]
  print(paste("Calculated ANOVA and saved plot for", expo))
}





## Calculate correlations (more carefully)
AC_labels <- (AC_PXS_means %>% drop_na())$assessment_center
AC_PXS_means <- AC_PXS_means %>% drop_na() %>% select(-assessment_center)
PXS_labels <- ukb_dict[match(pheno_list,ukb_dict$field),]$fieldname
AC_expo_means <- AC_expo_means %>% drop_na() %>% select(-assessment_center)
expos_labels <- ukb_dict[match(exposures_list,ukb_dict$field),]$fieldname

AC_PXS_corrs <- tibble(
  field1 = as.character(),
  field2 = as.character(),
  r = as.numeric(),
  p = as.numeric()
)
AC_expo_corrs <- AC_PXS_corrs

AC_PXS_corrs <- calculate_correlations(AC_PXS_means, AC_PXS_corrs)
AC_expo_corrs <- calculate_correlations(AC_expo_means, AC_expo_corrs)


n <- nrow(AC_PXS_corrs %>% filter(p_adj < 0.05))
N <- nrow(AC_PXS_corrs)
print(paste("For PXS, there are",n,"out of",N,"significantly correlated unique pairs"))

n <- nrow(AC_expo_corrs %>% filter(p_adj < 0.05))
N <- nrow(AC_expo_corrs)
print(paste("For exposures, there are",n,"out of",N,"significantly correlated unique pairs"))

# # PXSs PCA
# matrix <- AC_PXS_means
# colnames(matrix) <- PXS_labels
# matrix.pca <- prcomp(matrix, center=TRUE, scale. = TRUE)
# 
# matrix_t <- as_tibble(t(matrix))
# colnames(matrix_t) <- AC_labels
# matrix_t.pca <- prcomp(matrix_t, center=TRUE, scale. = TRUE)
# custom_ggbiplot(matrix.pca, labels = AC_labels, var.axes=FALSE)
# custom_ggbiplot(matrix_t.pca, labels = PXS_labels, var.axes=FALSE,
#                 labels.size = 3)
# 
# # exposures PXS
# matrix <- AC_expo_means
# colnames(matrix) <- expos_labels
# matrix.pca <- prcomp(matrix, center=TRUE, scale. = TRUE)
# 
# matrix_t <- as_tibble(t(matrix))
# colnames(matrix_t) <- AC_labels
# matrix_t.pca <- prcomp(matrix_t, center=TRUE, scale. = TRUE)
# custom_ggbiplot(matrix.pca, labels = AC_labels, var.axes=FALSE)
# custom_ggbiplot(matrix_t.pca, labels = expos_labels, var.axes=FALSE,
#                 labels.size = 2)

ggplot(AC_expo_means, aes(x=mean_f.1438.0.0, y=mean_f.1488.0.0)) +
  geom_text(aes(label=AC_labels))
ggplot(AC_PXS_means, aes(x=mean_PXS_AF, y=mean_PXS_COPD)) +
  geom_text(aes(label=AC_labels))

## Append field descriptions
ANOVA_PXS_tbl$field <- unname(sapply(ANOVA_PXS_tbl$field, function(x) str_split(x,"PXS_")[[1]][length(str_split(x,"PXS_")[[1]])], simplify = TRUE))
ANOVA_PXS_tbl <- ANOVA_PXS_tbl %>% left_join(ukb_dict, by="field")
ANOVA_expo_tbl <- ANOVA_expo_tbl %>% left_join(ukb_dict, by="field")
