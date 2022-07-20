## Libraries and directories ##
library(tidyverse)
library(data.table)
library(GGally)

dir_scratch <- "~/scratch3/PXS_pipeline/"
dir_script <- "~/jobs/PXS_pipeline/"
dir_data_showcase <- "~/scratch3/key_data/"

## Functions ##
calculate_ANOVA <- function(field, loc_out, ANOVA_tbl, AC_means, remove_negs=FALSE) {
  
  ANOVA_tibble <- pheno %>% select(assessment_center,all_of(field)) %>%
    rename(field_value = all_of(field)) %>%
    drop_na()
  if (remove_negs) {ANOVA_tibble <- ANOVA_tibble %>% filter(field_value >= 0)}
  # normalizes data since that is one of the ANOVA assumptions
  ANOVA_tibble$field_value <- qnorm(rank(ANOVA_tibble$field_value) / length(ANOVA_tibble$field_value))
  ANOVA_tibble <- ANOVA_tibble %>% filter(!is.infinite(field_value))
  
  aov1 <- summary(aov(data=ANOVA_tibble, field_value ~ ANOVA_tibble[["assessment_center"]]))
  f_stat <- aov1[[1]]$`F value`[1]
  p_value <- aov1[[1]]$`Pr(>F)`[1]
  
  
  ANOVA_tibble_summary <- ANOVA_tibble %>%
    drop_na() %>%
    group_by(assessment_center) %>%
    summarise(mean_field_value = mean(field_value, na.rm=TRUE),
              med_field_value = median(field_value, na.rm=TRUE),
              n = n()) %>%
    arrange(-mean_field_value)
  
  ANOVA_tbl <- ANOVA_tbl %>%
    add_row(
      field = field,
      f_stat = f_stat,
      p_value = p_value
    )
  
  # gg<-ggplot(ANOVA_tibble, aes(x=assessment_center,y=field_value)) +
  #   geom_boxplot() +
  #   geom_label(data=ANOVA_tibble_summary, aes(y=med_field_value, label=round(med_field_value,2)),
  #              label.padding = unit(0.1, "lines")) +
  #   coord_flip() +
  #   xlab("Assessment Center") +
  #   ylab(field) +
  #   labs(title = paste("Boxplots of",field,"per Assessment Center"),
  #        subtitle = paste0("ANOVA: F = ",round(f_stat,3),", p-value = ",round(p_value,3)))
  # ggsave(loc_out, gg, width = 3600, height = 2700, units = "px")
  
  ANOVA_tibble_summary <- ANOVA_tibble_summary %>%
    select(assessment_center, mean_field_value)
  colnames(ANOVA_tibble_summary)[2] <- paste0("mean_",field)
  AC_means <- AC_means %>% left_join(ANOVA_tibble_summary, by="assessment_center")
  
  out_list <- list(ANOVA_tbl, AC_means)
  out_list
}
calculate_correlations <- function(AC_means, AC_corrs) {
  for (i in 1:(length(colnames(AC_means))-1)) {
    col1 <- colnames(AC_means)[i]
    field1 <- substring(col1,6,nchar(col1))
    for (j in (i+1):length(colnames(AC_means)) ) {
      col2 <- colnames(AC_means)[j]
      field2 <- substring(col2,6,nchar(col2))
      
      cor1 <- cor.test(AC_means[[col1]], AC_means[[col2]])
      r = cor1$estimate[[1]]
      p = cor1$p.value
      
      AC_corrs <- AC_corrs %>%
        add_row(
          field1 = field1,
          field2 = field2,
          r = r,
          p = p
        )
      print(paste("Calculated correlation for", field1, "and", field2))
    }
  }
  AC_corrs$field1 <- unname(sapply(AC_corrs$field1, function(x) str_split(x,"PXS_")[[1]][length(str_split(x,"PXS_")[[1]])], simplify = TRUE))
  AC_corrs$field2 <- unname(sapply(AC_corrs$field2, function(x) str_split(x,"PXS_")[[1]][length(str_split(x,"PXS_")[[1]])], simplify = TRUE))
  AC_corrs <- AC_corrs %>%
    mutate(p_adj = p.adjust(AC_corrs$p, method="fdr")) %>%
    left_join(ukb_dict, by=c("field1"="field")) %>%
    rename(fieldname1 = fieldname) %>%
    left_join(ukb_dict, by=c("field2"="field")) %>%
    rename(fieldname2 = fieldname) %>%
    arrange(p_adj)
  AC_corrs
}






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
  rename(sex = f.31.0.0,
         year_of_birth = f.34.0.0,
         assessment_center = f.54.0.0) %>%
  mutate(sex = as.factor(sex),
         assessment_center = factor(assessment_center,
                                    levels = AC_tbl$Value,
                                    labels = AC_tbl$Meaning))

loc_phenolist <- paste0(dir_script,"phenotypes_ALL.txt")
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

ggpairs(AC_PXS_means)



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
AC_PXS_means <- AC_PXS_means %>% drop_na() %>% select(-assessment_center)
AC_expo_means <- AC_expo_means %>% drop_na() %>% select(-assessment_center)

AC_PXS_corrs <- tibble(
  field1 = as.character(),
  field2 = as.character(),
  r = as.numeric(),
  p = as.numeric()
)
AC_expo_corrs <- AC_PXS_corrs

AC_PXS_corrs <- calculate_correlations(AC_PXS_means, AC_PXS_corrs)
AC_expo_corrs <- calculate_correlations(AC_expo_means, AC_expo_corrs)


## Append field descriptions
ANOVA_PXS_tbl <- ANOVA_PXS_tbl %>% left_join(ukb_dict, by="field")
ANOVA_expo_tbl <- ANOVA_expo_tbl %>% left_join(ukb_dict, by="field")
