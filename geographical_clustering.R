## Libraries and directories ##
library(tidyverse)
library(data.table)

dir_scratch <- "~/scratch3/PXS_pipeline/"
dir_script <- "~/jobs/PXS_pipeline/"
dir_data_showcase <- "~/scratch3/key_data/"

## Code ##
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

#loc_fields <- paste0(dir_scratch,"fields_tbl.txt")
#fields <- as_tibble(fread(loc_fields))


loc_phenolist <- paste0(dir_script,"phenotypes_ALL.txt")
pheno_list <- readLines(loc_phenolist)

# Loops through fields and calculates divergence by
ANOVA_PXS_tbl <- tibble(
  field = as.character(),
  f_stat = as.numeric(),
  p_value = as.numeric()
)
for (field in pheno_list) {
  
  col_PXS <- paste0("PXS_",field)
  ANOVA_tibble <- pheno %>% select(assessment_center,all_of(col_PXS)) %>%
    rename(PXS = all_of(col_PXS)) %>%
    drop_na()
  
  aov1 <- summary(aov(data=ANOVA_tibble, PXS ~ ANOVA_tibble[["assessment_center"]]))
  f_stat <- aov1[[1]]$`F value`[1]
  p_value <- aov1[[1]]$`Pr(>F)`[1]
  
  
  ANOVA_tibble_summary <- ANOVA_tibble %>%
    drop_na() %>%
    group_by(assessment_center) %>%
    summarise(mean_PXS = mean(PXS, na.rm=TRUE),
              med_PXS = median(PXS, na.rm=TRUE),
              n = n()) %>%
    arrange(-mean_PXS)
  
  gg<-ggplot(ANOVA_tibble, aes(x=assessment_center,y=PXS)) +
    geom_boxplot() +
    geom_label(data=ANOVA_tibble_summary, aes(y=med_PXS, label=round(med_PXS,2)),
               label.padding = unit(0.1, "lines")) +
    coord_flip() +
    xlab("Assessment Center") +
    ylab(paste("PXS for",field)) +
    labs(title = paste("Boxplots of",field,"PXS per Assessment Center"),
         subtitle = paste0("ANOVA: F = ",round(f_stat,3),", p-value = ",round(p_value,3)))
  
  loc_out <- paste0(dir_scratch,field,"/","ANOVAbp_",field,".png")
  ggsave(loc_out, gg, width = 3600, height = 2700, units = "px")

  ANOVA_PXS_tbl <- ANOVA_PXS_tbl %>%
    add_row(
      field = field,
      f_stat = f_stat,
      p_value = p_value
    )
  print(paste("Calculated ANOVA and saved plot for", field))
}
ukb_dict <- as_tibble(fread(paste0(dir_data_showcase,"Data_Dictionary_Showcase.tsv"))) %>%
  mutate(field = paste0("f.",FieldID,".0.0")) %>%
  select(field, fieldname=Field) %>%
  add_row(
    field = c("AF", "CAD", "COPD", "T2D","fev1_inst1"),
    fieldname = c("Atrial fibrillation", "Coronary artery disease", "Chronic obstructive pulmonary disease",
                  "Type 2 diabetes", "Forced expiratory volume in 1-second (FEV1)")
  )
ANOVA_PXS_tbl <- ANOVA_PXS_tbl %>% left_join(ukb_dict, by="field")
