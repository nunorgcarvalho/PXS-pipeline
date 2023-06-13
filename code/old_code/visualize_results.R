# Need at least 8GB of RAM to load full LMM data later

# reuses code from aggregate_results and geographical_clustering to create clean
# tables and figures for publication

#### Libraries and directories ####
library(tidyverse)
library(data.table)
library(ggpubr)
library(ggrepel)
library(webshot)
library(gt)
# Functions #
source("~/jobs/PXS_pipeline/code/helper_functions.R")

dir_script <- "~/jobs/PXS_pipeline/code/"
dir_scratch <- "~/scratch3/PXS_pipeline/"
#dir_scratch <- "~/scratch3/08-01_PXS_pipeline/"
dir_data_showcase <- "~/scratch3/key_data/" # contains 'Data_Dictionary_Showcase.tsv' from UKBB

setwd(dir_scratch)

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

mean(REML_expo$h2g) + 1.96*sqrt(sum(REML_expo$h2g_err**2))


table1 <- left_join(REML_expo, LMM_expo, by="field") %>%
  select(fieldname, h2g, h2g_err, lambda, N_above_bonferroni) %>%
  mutate(h2g_CI = paste0("(",round(h2g - 1.96*h2g_err,rounding_decimals),", ",
                         round(h2g + 1.96*h2g_err,rounding_decimals),")")) %>%
  arrange(-h2g) %>%
  drop_na() %>%
  mutate(n = row_number()) %>%
  select(n,fieldname, h2g, h2g_CI, lambda, N_above_bonferroni) %>%
  filter(row_number() <= 20)

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
    footnote = "P < 0.05 / 19400443 (number of total SNPs); lambda-adjusted",
    locations = cells_column_labels(columns = c("N_above_bonferroni"))
  )
  
table1gt
gtsave(table1gt, "table1_20.png","./figures/")

#### Table 2 ####
REML_PXS <- as_tibble(fread(paste0(dir_scratch,"REML_PXS_results.txt")))
LMM_PXS <- as_tibble(fread(paste0(dir_scratch,"LMM_PXS_results.txt"))) %>%
  select(-fieldname)


table2 <- left_join(REML_PXS, LMM_PXS, by="field") %>%
  select(fieldname, h2g, h2g_err, lambda, N_above_bonferroni) %>%
  mutate(h2g_CI = paste0("(",round(h2g - 1.96*h2g_err,rounding_decimals),", ",
                         round(h2g + 1.96*h2g_err,rounding_decimals),")")) %>%
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
    title = "Heritability and Lambda Inflation Factor of PXSs",
    subtitle = paste0("Displaying ", nrow(table2)," cardiometabolic health outcome PXSs")
  ) %>%
  tab_footnote(
    footnote = "P < 0.05 / 19400443 (number of total SNPs); lambda-adjusted",
    locations = cells_column_labels(columns = c("N_above_bonferroni"))
  )

table2gt
gtsave(table2gt, "table2.png","./figures/")

#### Table 3 ####
genCorrs <- as_tibble(fread(paste0(dir_scratch,"genCorr_CRF_results.txt")))


table3 <- genCorrs %>%
  select(CRF_fieldname, h2g2, h2g2_err, gencorr, gencorr_err) %>%
  mutate(h2g2_CI = paste0("(",round(h2g2 - 1.96*h2g2_err,rounding_decimals),", ",
                         round(h2g2 + 1.96*h2g2_err,rounding_decimals),")"),
         gencorr_CI = paste0("(",round(gencorr - 1.96*gencorr_err,rounding_decimals),", ",
                          round(gencorr + 1.96*gencorr_err,rounding_decimals),")")) %>%
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
    subtitle = paste0("Displaying ",nrow(table3)," Clinical Risk Factors of Type 2 Diabetes")
  )

table3gt
gtsave(table3gt, "table3.png","./figures/")

# Manhattan Plot settings
cols_to_keep <- c("SNP","CHR","BP","BETA","P_BOLT_LMM_INF","CHISQ_BOLT_LMM_INF")

width = 7200
height = 3200

#### Figure 2 ####
expos_to_plot <- c("f.20116.0.0","f.6139.0.0")
plots2 <- list()
for (i in 1:length(expos_to_plot)) {
  expo <- expos_to_plot[i]
  fieldname <- (REML_expo %>% filter(field==expo))$fieldname[1]
  loc_LMM <- paste0(dir_scratch,"exposures/",expo,"/LMM_",expo,"_bgen.txt")
  if (!file.exists(loc_LMM)) {next}
  LMM <- as_tibble(fread(loc_LMM, select=cols_to_keep)) %>%
    rename(P_unadj = P_BOLT_LMM_INF)
  
  N <- nrow(LMM)
  bonferroni <- 0.05 / N
  lambda <- get_lambda(LMM$CHISQ_BOLT_LMM_INF)
  h2 <- (REML_expo %>% filter(field==expo))$h2g[[1]]
  
  # adjusts for lambda inflation factor
  LMM <- LMM %>% mutate(
    P = exp(pchisq(CHISQ_BOLT_LMM_INF/lambda,1, lower.tail=FALSE, log.p=TRUE))
  )
  
  # reduces the number of points to plot significantly while trying to keep
  # the Manhattan plot visually the same as before
  LMM2 <- downscale_sf(LMM)
  
  gg <- plot_Manhattan(LMM2, fieldname, N, lambda, h2)
  plots2[[i]] <- gg
  
  print(paste("Saved Manhattan plot for",expo))
}

fig2 <- ggarrange(plots2[[1]],plots2[[2]])
loc_out <- paste0("./figures/","figure2.png")
ggsave(loc_out,plot=fig2,width = width, height = height, units="px")

#### Figure 3 ####
phenos_to_plot <- c("AF","f.30740.0.0")
plots3 <- list()
for (i in 1:length(phenos_to_plot)) {
  pheno <- phenos_to_plot[i]
  fieldname <- paste0((REML_PXS %>% filter(field==pheno))$fieldname[1], " PXS")
  loc_LMM <- paste0(dir_scratch,pheno,"/LMM_",pheno,"_bgen.txt")
  if (!file.exists(loc_LMM)) {next}
  LMM <- as_tibble(fread(loc_LMM, select=cols_to_keep)) %>%
    rename(P_unadj = P_BOLT_LMM_INF)
  
  N <- nrow(LMM)
  bonferroni <- 0.05 / N
  lambda <- get_lambda(LMM$CHISQ_BOLT_LMM_INF)
  h2 <- (REML_PXS %>% filter(field==pheno))$h2g[[1]]
  
  # adjusts for lambda inflation factor
  LMM <- LMM %>% mutate(
    P = exp(pchisq(CHISQ_BOLT_LMM_INF/lambda,1, lower.tail=FALSE, log.p=TRUE))
  )
  
  # reduces the number of points to plot significantly while trying to keep
  # the Manhattan plot visually the same as before
  LMM2 <- downscale_sf(LMM)
  
  gg <- plot_Manhattan(LMM2, fieldname, N, lambda, h2)
  plots3[[i]] <- gg
  
  print(paste("Saved Manhattan plot for",pheno))
}

fig3 <- ggarrange(plots3[[1]],plots3[[2]])
loc_out <- paste0("./figures/","figure3.png")
ggsave(loc_out,plot=fig3,width = width, height = height, units="px")


#### ANOVA-related settings and setup ####
loc_phenolist <- paste0(dir_script,"../input_data/phenotypes.txt")
pheno_list <- readLines(loc_phenolist)

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

ANOVA_expo_tbl <- tibble(
  field = as.character(),
  f_stat = as.numeric(),
  p_value = as.numeric()
)
AC_expo_means <- tibble(assessment_center = levels(pheno$assessment_center))

#### Table 4 ####

# Loops through exposures
loc_fields <- paste0(dir_scratch,"fields_tbl.txt")
fields <- as_tibble(fread(loc_fields))

# only performs analysis on quantitative exposures found in the data
exposures_list_quant <- (fields %>% filter(use_type == "exposure",
                                     data_type == "quantitative",
                                     ))$field
exposures_list_quant <- exposures_list_quant[exposures_list_quant %in% colnames(pheno)]

expo_to_plot <- "f.1070.0.0"

for (expo in exposures_list_quant) {
  dir_expo <-  paste0(dir_scratch,"exposures/",expo)
  loc_out <- paste0(dir_expo,"/","ANOVAbp_",expo,".png")
  out_list <- calculate_ANOVA(expo, loc_out, ANOVA_expo_tbl, AC_expo_means, TRUE, FALSE)
  ANOVA_expo_tbl <- out_list[[1]]
  AC_expo_means <- out_list[[2]]
  if (expo %in% expo_to_plot) {
    fig4tbl1 <- out_list[[3]]
    fig4tbl2 <- out_list[[4]]
  }
  print(paste("Calculated ANOVA for", expo))
}
ANOVA_expo_tbl <- ANOVA_expo_tbl %>% left_join(ukb_dict, by="field")

table4 <- ANOVA_expo_tbl %>%
  mutate(p_text = ifelse(p_value < 1e-320,"< 1E-320",paste0(formatC(p_value,digits=2,format="E")))) %>%
  arrange(-f_stat) %>%
  mutate(n = row_number()) %>%
  select(n,fieldname, f_stat, p_text) %>%
  filter(row_number() <= 10)

table4gt <- gt(table4) %>%
  gt_theme() %>%
  cols_label(
    n = "#",
    fieldname = "Exposure",
    f_stat = "F-statistic",
    p_text = "p-value"
  ) %>%
  fmt_number(
    columns = c("f_stat"),
    suffixing = F,
    decimals = rounding_decimals
  ) %>%
  cols_align(
    align = "left",
    columns = c("fieldname")
  ) %>%
  tab_header(
    title = "Geospatial Heterogeneity (ANOVA) of Exposures by Assessment Center",
    subtitle = paste0("Displaying top ",nrow(table4)," exposures by F-statistic")
  )

table4gt
gtsave(table4gt, "table4.png","./figures/")

#### Figure 4 ####

the_fieldname <- (REML_expo %>% filter(field==expo_to_plot))$fieldname[1]
f_stat <- (table4 %>% filter(fieldname==the_fieldname))$f_stat[1]
p_text <- (table4 %>% filter(fieldname==the_fieldname))$p_text[1]
fig4 <- ggplot(fig4tbl1, aes(x=assessment_center,y=field_value)) +
  geom_boxplot() +
  geom_label(data=fig4tbl2, aes(y=med_field_value, label=round(med_field_value,2)),
             label.padding = unit(0.1, "lines")) +
  coord_flip() +
  xlab("Assessment Center") +
  ylab(paste0("Phenotype Value (z-score normalized)")) +
  labs(title = paste("Boxplots of",the_fieldname,"per Assessment Center"),
       subtitle = paste0("ANOVA: F = ",round(f_stat,1)," :: p-value ",p_text)) +
  theme_light()
loc_out <- paste0("./figures/","figure4.png")
ggsave(loc_out,plot=fig4,width = width/2, height = height, units="px")

#### Table 6 ####

ANOVA_PXS_tbl <- tibble(
  field = as.character(),
  f_stat = as.numeric(),
  p_value = as.numeric()
)
AC_PXS_means <- tibble(assessment_center = levels(pheno$assessment_center))

PXS_to_plot = "CAD"
for (field in pheno_list) {
  
  col_PXS <- paste0("PXS_",field)
  loc_out <- paste0(dir_scratch,field,"/","ANOVAbp_",field,".png")
  out_list <- calculate_ANOVA(col_PXS, loc_out, ANOVA_PXS_tbl, AC_PXS_means, FALSE, FALSE)
  ANOVA_PXS_tbl <- out_list[[1]]
  AC_PXS_means <- out_list[[2]]
  if (field == PXS_to_plot) {
    fig6tbl1 <- out_list[[3]]
    fig6tbl2 <- out_list[[4]]
  }
  print(paste("Calculated ANOVA for", field))
}
ANOVA_PXS_tbl$field <- unname(sapply(ANOVA_PXS_tbl$field, function(x) str_split(x,"PXS_")[[1]][length(str_split(x,"PXS_")[[1]])], simplify = TRUE))
ANOVA_PXS_tbl <- ANOVA_PXS_tbl %>% left_join(ukb_dict, by="field")

table6 <- ANOVA_PXS_tbl %>%
  mutate(p_text = ifelse(p_value < 1e-320,"< 1E-320",paste0(formatC(p_value,digits=2,format="E")))) %>%
  arrange(-f_stat) %>%
  mutate(n = row_number()) %>%
  select(n,fieldname, f_stat, p_text)

table6gt <- gt(table6) %>%
  gt_theme() %>%
  cols_label(
    n = "#",
    fieldname = "PXS",
    f_stat = "F-statistic",
    p_text = "p-value"
  ) %>%
  fmt_number(
    columns = c("f_stat"),
    suffixing = F,
    decimals = rounding_decimals
  ) %>%
  cols_align(
    align = "left",
    columns = c("fieldname")
  ) %>%
  tab_header(
    title = "Geospatial Heterogeneity (ANOVA) of PXSs by Assessment Center",
    subtitle = paste0("Displaying ",nrow(table6)," cardiometabolic health outcomes PXSs")
  )

table6gt
gtsave(table6gt, "table6.png","./figures/")

#### Figure 6 ####

the_fieldname <- (REML_PXS %>% filter(field==PXS_to_plot))$fieldname[1]
f_stat <- (table6 %>% filter(fieldname==the_fieldname))$f_stat[1]
p_text <- (table6 %>% filter(fieldname==the_fieldname))$p_text[1]
fig6 <- ggplot(fig6tbl1, aes(x=assessment_center,y=field_value)) +
  geom_boxplot() +
  geom_label(data=fig6tbl2, aes(y=med_field_value, label=round(med_field_value,2)),
             label.padding = unit(0.1, "lines")) +
  coord_flip() +
  xlab("Assessment Center") +
  ylab(paste0("Phenotype Value (z-score normalized)")) +
  labs(title = paste("Boxplots of",the_fieldname,"PXS per Assessment Center"),
       subtitle = paste0("ANOVA: F = ",round(f_stat,1)," :: p-value ",p_text)) +
  theme_light()
loc_out <- paste0("./figures/","figure6.png")
ggsave(loc_out,plot=fig6,width = width/2, height = height, units="px")

AC_labels <- (AC_PXS_means %>% drop_na())$assessment_center
#### Table 5 ####
AC_expo_means <- AC_expo_means %>% drop_na() %>% select(-assessment_center)
expos_labels <- ukb_dict[match(exposures_list_quant,ukb_dict$field),]$fieldname

AC_expo_corrs <- tibble(
  field1 = as.character(),
  field2 = as.character(),
  r = as.numeric(),
  p = as.numeric()
)
AC_expo_corrs <- calculate_correlations(AC_expo_means, AC_expo_corrs)
n <- nrow(AC_expo_corrs %>% filter(p_adj < 0.05))
N <- nrow(AC_expo_corrs)
print(paste("For exposures, there are",n,"out of",N,"significantly correlated unique pairs"))

table5 <- AC_expo_corrs %>%
  arrange(p_adj) %>%
  mutate(n = row_number(),
         p_text = ifelse(p_adj < 1e-320,"< 1E-320",paste0(formatC(p_adj,digits=2,format="E")))) %>%
  select(n, fieldname1, fieldname2, r, p_text) %>%
  filter(row_number() <= 10)

table5gt <- gt(table5) %>%
  gt_theme() %>%
  cols_label(
    n = "#",
    fieldname1 = "Exposure 1",
    fieldname2 = "Exposure 2",
    r = "correlation",
    p_text = "p-value"
  ) %>%
  fmt_number(
    columns = c("r"),
    suffixing = F,
    decimals = rounding_decimals
  ) %>%
  cols_align(
    align = "left",
    columns = c("fieldname1","fieldname2")
  ) %>%
  tab_header(
    title = "Geospatial Correlations of Exposures across Assessment Center",
    subtitle = paste0("Displaying top ",nrow(table5)," exposure pairs by p-value. ",
                      n," out of ",N," (",round(n/N,3)*100,"%) unique exposure pairs are significantly correlated")
  ) %>%
  tab_footnote(
    footnote = "FDR-adjusted",
    locations = cells_column_labels(columns = c("p_text"))
  )

table5gt
gtsave(table5gt, "table5.png","./figures/")

#### Table 7 ####
AC_PXS_means <- AC_PXS_means %>% drop_na() %>% select(-assessment_center)
PXS_labels <- ukb_dict[match(pheno_list,ukb_dict$field),]$fieldname

AC_PXS_corrs <- tibble(
  field1 = as.character(),
  field2 = as.character(),
  r = as.numeric(),
  p = as.numeric()
)

AC_PXS_corrs <- calculate_correlations(AC_PXS_means, AC_PXS_corrs)

n <- nrow(AC_PXS_corrs %>% filter(p_adj < 0.05))
N <- nrow(AC_PXS_corrs)
print(paste("For PXS, there are",n,"out of",N,"significantly correlated unique pairs"))

table7 <- AC_PXS_corrs %>%
  arrange(p_adj) %>%
  mutate(n = row_number(),
         p_text = ifelse(p_adj < 1e-320,"< 1E-320",paste0(formatC(p_adj,digits=2,format="E")))) %>%
  select(n, fieldname1, fieldname2, r, p_text) %>%
  filter(row_number() <= 10)

table7gt <- gt(table7) %>%
  gt_theme() %>%
  cols_label(
    n = "#",
    fieldname1 = "PXS 1",
    fieldname2 = "PXS 2",
    r = "correlation",
    p_text = "p-value"
  ) %>%
  fmt_number(
    columns = c("r"),
    suffixing = F,
    decimals = rounding_decimals
  ) %>%
  cols_align(
    align = "left",
    columns = c("fieldname1","fieldname2")
  ) %>%
  tab_header(
    title = "Geospatial Correlations of PXSs across Assessment Center",
    subtitle = paste0("Displaying top ",nrow(table5)," PXS pairs by p-value. ",
                      n," out of ",N," (",round(n/N,3)*100,"%) unique PXS pairs are significantly correlated")
  ) %>%
  tab_footnote(
    footnote = "FDR-adjusted",
    locations = cells_column_labels(columns = c("p_text"))
  )

table7gt
gtsave(table7gt, "table7.png","./figures/")

#### Figure 5 ####

pairs_to_plot <- list(c("f.1050.0.0","f.1070.0.0"),c("f.21001.0.0","f.30760.0.0"))

fig5tbl <- AC_expo_means %>%
  rename(x = all_of(paste0("mean_",pairs_to_plot[[1]][1])), y = all_of(paste0("mean_",pairs_to_plot[[1]][2]))) %>%
  mutate(assessment_center = AC_labels) %>%
  select(assessment_center, x, y)
fig5 <- ggplot(fig5tbl, aes(x, y)) +
  geom_smooth(method = "lm", color="#5BC6DD") +
  geom_point() +
  geom_text_repel(aes(label=assessment_center), seed=1) +
  xlab(paste0("Mean ",(REML_expo %>% filter(field==pairs_to_plot[[1]][1]))$fieldname[1], " (z-score normalized)")) +
  ylab(paste0("Mean ",(REML_expo %>% filter(field==pairs_to_plot[[1]][2]))$fieldname[1], " (z-score normalized)")) +
  labs(title = "Geospatial Correlations of Two Exposures across Assessment Center",
       subtitle = paste0("r = ", round((AC_expo_corrs %>% filter(field1==pairs_to_plot[[1]][1],field2==pairs_to_plot[[1]][2]))$r[1], rounding_decimals),
                         " :: p-value = ", formatC((AC_expo_corrs %>% filter(field1==pairs_to_plot[[1]][1],field2==pairs_to_plot[[1]][2]))$p_adj[1],format="E",digits=rounding_decimals))) +
  theme_light()

loc_out <- paste0("./figures/","figure5.png")
ggsave(loc_out,plot=fig5,width = (width/2), height = height, units="px")

#### Figure 7 ####

fig7tbl <- AC_PXS_means %>%
  rename(x = all_of(paste0("mean_PXS_",pairs_to_plot[[2]][1])), y = all_of(paste0("mean_PXS_",pairs_to_plot[[2]][2]))) %>%
  mutate(assessment_center = AC_labels) %>%
  select(assessment_center, x, y)
fig7 <- ggplot(fig7tbl, aes(x, y)) +
  geom_smooth(method = "lm", color="#5BC6DD") +
  geom_point() +
  geom_text_repel(aes(label=assessment_center), seed=1) +
  xlab(paste0("Mean ",(REML_PXS %>% filter(field==pairs_to_plot[[2]][1]))$fieldname[1], " (z-score normalized)")) +
  ylab(paste0("Mean ",(REML_PXS %>% filter(field==pairs_to_plot[[2]][2]))$fieldname[1], " (z-score normalized)")) +
  labs(title = "Geospatial Correlations of Two PXSs across Assessment Center",
       subtitle = paste0("r = ", round((AC_PXS_corrs %>% filter(field1==pairs_to_plot[[2]][1],field2==pairs_to_plot[[2]][2]))$r[1], rounding_decimals),
                         " :: p-value = ", formatC((AC_PXS_corrs %>% filter(field1==pairs_to_plot[[2]][1],field2==pairs_to_plot[[2]][2]))$p_adj[1],format="E",digits=rounding_decimals))) +
  theme_light()

loc_out <- paste0("./figures/","figure7.png")
ggsave(loc_out,plot=fig7,width = (width/2), height = height, units="px")

#### Table 8 ####
envLM <- as_tibble(fread(paste0(dir_scratch,"envLM_PXS_results.txt")))

table8 <- envLM %>%
  arrange(-R2adj) %>%
  mutate(n = row_number(),
         R2_text = paste0(round(R2adj*100,2),"%")) %>%
  select(n,fieldname, R2_text)

table8gt <- gt(table8) %>%
  gt_theme() %>%
  cols_label(
    n = "#",
    fieldname = "PXS",
    R2_text = md("R<sup>2</sup> (adj)")
  ) %>%
  cols_align(
    align = "left",
    columns = c("fieldname")
  ) %>%
  tab_header(
    title = "Environmental Linear Model for PXS",
    subtitle = paste0("PXS ~ sex + age + assessment_center + 40 PCs")
  ) %>%
  tab_footnote(
    footnote = "All are significant with P < 1E-320",
    locations = cells_column_labels(columns = c("R2_text"))
  )

table8gt
gtsave(table8gt, "table8.png","./figures/")
