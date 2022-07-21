# Need at least 8GB of RAM to load full LMM data later
## Libraries and directories ##
library(tidyverse)
library(data.table)
library(ggrepel)

dir_scratch <- "~/scratch3/PXS_pipeline/"
dir_script <- "~/jobs/PXS_pipeline/"
dir_data_showcase <- "~/scratch3/key_data/"

## Functions ##
extract_from_REML <- function(REML) {
  N <- length(REML)
  
  # extract time for analysis
  line <- REML[N]
  elapsed_hours <- as.numeric(str_split(line," ")[[1]][7]) / (60 * 60)
  
  # extract h2e
  line <- REML[N-4]
  line_split <- str_split(line, ": ")[[1]][2]
  h2e <- as.numeric(str_split(line_split, " ")[[1]][1])
  h2e_err <- parse_number(str_split(line_split, " ")[[1]][2])
  
  # extract h2g
  line <- REML[N-2]
  line_split <- str_split(line, ": ")[[1]][2]
  h2g <- as.numeric(str_split(line_split, " ")[[1]][1])
  h2g_err <- parse_number(str_split(line_split, " ")[[1]][2])
  
  out <- c(h2e, h2e_err, h2g, h2g_err, elapsed_hours)
  out
}

## Code ##

### REML results for each phenotypes
loc_phenolist <- paste0(dir_script,"phenotypes_ALL.txt")
pheno_list <- readLines(loc_phenolist)

REML_PXS_tbl <- tibble(
  field = as.character(),
  h2e = as.numeric(),
  h2e_err = as.numeric(),
  h2g = as.numeric(),
  h2g_err = as.numeric(),
  elapsed_hours = as.numeric()
)
REML_expo_tbl <- REML_PXS_tbl
for (pheno in pheno_list) {
  loc_REML <- paste0(dir_scratch,pheno,"/",pheno,"_PXS_BOLTREML.out")
  REML <- readLines(loc_REML)
  out <- extract_from_REML(REML)
  
  # adds row to table
  REML_PXS_tbl <- REML_PXS_tbl %>%
    add_row(
      field = pheno,
      h2e = out[1],
      h2e_err = out[2],
      h2g = out[3],
      h2g_err = out[4],
      elapsed_hours = out[5]
    )
  
  print(paste("Read REML results for",pheno))
}

### REML results for each exposure
loc_expolist <- paste0(dir_script,"exposures.txt")
exposures_list <- readLines(loc_expolist)

for (expo in exposures_list) {
  loc_REML <- paste0(dir_scratch,"exposures/",expo,"/",expo,"_BOLTREML.out")
  if (!file.exists(loc_REML)) {next}
  REML <- readLines(loc_REML)
  out <- extract_from_REML(REML)
  
  # adds row to table
  REML_expo_tbl <- REML_expo_tbl %>%
    add_row(
      field = expo,
      h2e = out[1],
      h2e_err = out[2],
      h2g = out[3],
      h2g_err = out[4],
      elapsed_hours = out[5]
    )
  
  print(paste("Read REML results for",expo))
}

### Appends fields name to REML tables
ukb_dict <- as_tibble(fread(paste0(dir_data_showcase,"Data_Dictionary_Showcase.tsv"))) %>%
  mutate(field = paste0("f.",FieldID,".0.0")) %>%
  select(field, fieldname=Field) %>%
  add_row(
    field = c("AF", "CAD", "COPD", "T2D","fev1_inst1"),
    fieldname = c("Atrial fibrillation", "Coronary artery disease", "Chronic obstructive pulmonary disease",
                  "Type 2 diabetes", "Forced expiratory volume in 1-second (FEV1)")
  )
REML_PXS_tbl <- REML_PXS_tbl %>% left_join(ukb_dict, by="field")
REML_expo_tbl <- REML_expo_tbl %>% left_join(ukb_dict, by="field")

loc_out <- paste0(dir_scratch, "REML_PXS_results.txt")
write.table(REML_PXS_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

loc_out <- paste0(dir_scratch, "REML_exposures_results.txt")
write.table(REML_expo_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)


#### GWAS ANALYSIS ########

# Functions #
plot_Manhattan <- function(LMM, fieldname, N=-1, lambda=-1) {
  ### MY OWN MANHATTAN FUNCTION WITH HELP FROM https://r-graph-gallery.com/101_Manhattan_plot.html ###
  
  if (N==-1) {N <- nrow(LMM)}
  bonferroni <- 0.05 / N
  if (lambda==-1) {lambda <- get_lambda(LMM$CHISQ_BOLT_LMM_INF)}
  
  manhattan_tbl <- LMM %>%
    group_by(CHR) %>%
    summarize(chr_length = max(BP)) %>%
    mutate(total_bp = cumsum(as.numeric(chr_length)) - chr_length) %>%
    select(-chr_length) %>%
    left_join(LMM2, ., by="CHR") %>%
    arrange(CHR, BP) %>%
    mutate(BP_cum = total_bp + BP)
  
  axes_tbl <- manhattan_tbl %>%
    group_by(CHR) %>%
    summarize( center_bp = mean(c(max(BP_cum),min(BP_cum))) )
  
  if (sig_SNPs)
  sig_SNPs <- manhattan_tbl %>% filter(P < bonferroni) %>%
    group_by(CHR) %>% filter(P == min(P))
  max_logP <- ceiling(max(-log10(manhattan_tbl$P)))
  
  gg<-ggplot(manhattan_tbl, aes(x = BP_cum, y = -log10(P))) +
    geom_abline(slope=0,intercept = -log10(bonferroni), color = "red") +
    geom_point(aes(color=as.factor(CHR)), alpha = 0.8, size = 1.3) +
    geom_text_repel(data=sig_SNPs, aes(label=SNP,x=BP_cum,y=-log10(P)), nudge_y = 0.15, seed=1) +
    scale_color_manual(values = rep(c("skyblue","dodgerblue1"), 22)) +
    scale_x_continuous(label = axes_tbl$CHR, breaks = axes_tbl$center_bp) +
    scale_y_continuous(expand = c(0,0), breaks = c(0:max_logP), limits = c(0,max_logP+0.5)) +
    theme_light() +
    xlab("Chromosome") +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(title = paste0("Manhattan plot for ", fieldname),
         subtitle = paste0("Number of SNPs = ", round(N / 1E6,1), "M :: ",
                           "Number of significant SNPs = ", nrow(manhattan_tbl %>% filter(P<bonferroni)),
                           " :: Lambda = ", round(lambda,3)))
  gg
}
downscale_sf <- function(sf, base = 4) {
  for (i in c(-1,-2,-3)) {
    top <- 10**(i+1)
    bottom <- 10**(i)
    downscale <- base**(4+i)
    sf <- sf %>% filter( (P > top)|(P < bottom)|(row_number() %% downscale == 0) )
  }
  sf
}
get_lambda <- function(LMM_chisq) {
  lambda <- median(LMM_chisq) / qchisq(0.5,1)
  lambda
}

# Code
cols_to_keep <- c("SNP","CHR","BP","P_BOLT_LMM_INF","CHISQ_BOLT_LMM_INF")


LMM_PXS_tbl <- tibble(
  field = as.character(),
  lambda = as.numeric(),
  N = as.numeric(),
  N_above_bonferroni = as.numeric(),
  N_above_5E_8 = as.numeric()
)
LMM_expo_tbl <- LMM_PXS_tbl

# Loop through disease PXSs
for (pheno in pheno_list) {
  
  fieldname <- paste0(pheno, " PXS")
  loc_LMM <- paste0(dir_scratch,pheno,"/LMM_",pheno,"_bgen.txt")
  LMM <- as_tibble(fread(loc_LMM, select=cols_to_keep)) %>%
    rename(P = P_BOLT_LMM_INF)
  
  N <- nrow(LMM)
  lambda <- get_lambda(LMM$CHISQ_BOLT_LMM_INF)
  LMM_PXS_tbl <- LMM_PXS_tbl %>%
    add_row(
      field = pheno,
      lambda = lambda,
      N = N,
      N_above_bonferroni = nrow(LMM %>% filter(P < bonferroni)),
      N_above_5E_8 = nrow(LMM %>% filter(P < 5E-8))
    )
  
  
  # reduces the number of points to plot significantly while trying to keep
  # the Manhattan plot visually the same as before
  LMM2 <- downscale_sf(LMM)
  
  gg <- plot_Manhattan(LMM2, fieldname, N, lambda)
  
  loc_out <- paste0(dir_scratch,pheno,"/LMM_",pheno,"_manhattan.png")
  ggsave(loc_out, gg, width=3600,height=2700, units="px")
  
  print(paste("Saved Manhattan plot for",pheno))
}
# Loop through exposures
for (expo in exposures_list) {
  
  fieldname <- expo
  loc_LMM <- paste0(dir_scratch,"exposures/",expo,"/LMM_",expo,"_bgen.txt")
  if (!file.exists(loc_LMM)) {next}
  LMM <- as_tibble(fread(loc_LMM, select=cols_to_keep)) %>%
    rename(P = P_BOLT_LMM_INF)
  
  N <- nrow(LMM)
  lambda <- get_lambda(LMM$CHISQ_BOLT_LMM_INF)
  LMM_expo_tbl <- LMM_expo_tbl %>%
    add_row(
      field = expo,
      lambda = lambda,
      N = N,
      N_above_bonferroni = nrow(LMM %>% filter(P < bonferroni)),
      N_above_5E_8 = nrow(LMM %>% filter(P < 5E-8))
    )
  
  
  # reduces the number of points to plot significantly while trying to keep
  # the Manhattan plot visually the same as before
  LMM2 <- downscale_sf(LMM)
  
  gg <- plot_Manhattan(LMM2, fieldname, N, lambda)
  
  loc_out <- paste0(dir_scratch,"exposures/",expo,"/LMM_",expo,"_manhattan.png")
  ggsave(loc_out, gg, width=3600,height=2700, units="px")
  
  print(paste("Saved Manhattan plot for",expo))
}


LMM_PXS_tbl <- LMM_PXS_tbl %>% left_join(ukb_dict, by="field")
LMM_expo_tbl <- LMM_expo_tbl %>% left_join(ukb_dict, by="field")

loc_out <- paste0(dir_scratch, "LMM_PXS_results.txt")
write.table(LMM_PXS_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

loc_out <- paste0(dir_scratch, "LMM_exposures_results.txt")
write.table(LMM_expo_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)
