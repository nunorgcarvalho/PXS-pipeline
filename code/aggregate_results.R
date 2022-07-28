# Need at least 8GB of RAM to load full LMM data later
## Libraries and directories ##
library(tidyverse)
library(data.table)
library(ggrepel)

dir_script <- "~/jobs/PXS_pipeline/code/"
dir_scratch <- "~/scratch3/PXS_pipeline/"
dir_data_showcase <- "~/scratch3/key_data/" # contains 'Data_Dictionary_Showcase.tsv' from UKBB

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
extract_from_genCorr <- function(genCorr) {
  N <- length(genCorr)
  
  # extract time for analysis
  line <- genCorr[N]
  elapsed_hours <- as.numeric(str_split(line," ")[[1]][7]) / (60 * 60)
  
  # extract h2g1
  line <- genCorr[N-4]
  line_split <- str_split(line, ": ")[[1]][2]
  h2g1 <- as.numeric(str_split(line_split, " ")[[1]][1])
  h2g1_err <- parse_number(str_split(line_split, " ")[[1]][2])
  
  # extract gencorr
  line <- genCorr[N-3]
  line_split <- str_split(line, ": ")[[1]][2]
  gencorr <- as.numeric(str_split(line_split, " ")[[1]][1])
  gencorr_err <- parse_number(str_split(line_split, " ")[[1]][2])
  
  # extract h2g2
  line <- genCorr[N-2]
  line_split <- str_split(line, ": ")[[1]][2]
  h2g2 <- as.numeric(str_split(line_split, " ")[[1]][1])
  h2g2_err <- parse_number(str_split(line_split, " ")[[1]][2])
  
  out <- c(h2g1, h2g1_err, gencorr, gencorr_err, h2g2, h2g2_err, elapsed_hours)
  out
}
extract_from_envLM <- function(envLM) {
  N <- length(envLM)
  
  # extract R2
  line <- envLM[N-3]
  line_split <- str_split(substr(line,22,nchar(line)),",\tAdjusted R-squared:  ")
  R2 <- as.numeric(line_split[[1]][1])
  R2adj <- as.numeric(line_split[[1]][2])
  
  # extract F and p
  line <- envLM[N-2]
  line_split <- str_split(substr(line,14,nchar(line))," ")
  F_stat <- as.numeric(line_split[[1]][length(line_split[[1]]) - 9])
  p <- as.numeric(line_split[[1]][length(line_split[[1]])])
  
  out <- c(R2, R2adj, F_stat, p)
  out
}

## Code ##

### REML+envLM results for each phenotypes
loc_phenolist <- paste0(dir_script,"../input_data/phenotypes.txt")
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
envLM_PXS_tbl <- tibble(
  field = as.character(),
  R2 = as.numeric(),
  R2adj = as.numeric(),
  F_stat = as.numeric(),
  p = as.numeric()
)
for (pheno in pheno_list) {
  # REML
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
  
  # enVLM
  loc_envLM <- paste0(dir_scratch,pheno,"/",pheno,"_compute_PXS_LM.Rout")
  envLM <- readLines(loc_envLM)
  out <- extract_from_envLM(envLM)
  
  envLM_PXS_tbl <- envLM_PXS_tbl %>% add_row(
    field = pheno,
    R2 = out[1],
    R2adj = out[2],
    F_stat = out[3],
    p = out[4]
  )
  
  print(paste("Read envLM results for",pheno))
}

### REML results for each exposure
loc_expolist <- paste0(dir_script,"../input_data/exposures.txt")
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

### REML results for genCorr with CRFs
genCorr_CRF_tbl <- tibble(
  pheno_field = as.character(),
  h2g1 = as.numeric(),
  h2g1_err = as.numeric(),
  CRF_field = as.character(),
  h2g2 = as.numeric(),
  h2g2_err = as.numeric(),
  gencorr = as.numeric(),
  gencorr_err = as.numeric(),
  elapsed_hours = as.numeric()
)
for (pheno in pheno_list) {
  loc_CRFs <- paste0(dir_scratch,pheno,"/",pheno,"_CRFs.txt")
  if (!file.exists(loc_CRFs)) {next}
  CRFs_list <- readLines(loc_CRFs)
  for (CRF in CRFs_list) {
    #if (CRF == "f.4080.0.0") {next}
    loc_genCorr <- paste0(dir_scratch,pheno,"/PXS_",pheno,"_",CRF,"_genCorr.out")
    genCorr <- readLines(loc_genCorr)
    
    out <- extract_from_genCorr(genCorr)
    
    # adds row to table
    genCorr_CRF_tbl <- genCorr_CRF_tbl %>% add_row(
      pheno_field = pheno,
      h2g1 = out[1],
      h2g1_err = out[2],
      CRF_field = CRF,
      h2g2 = out[5],
      h2g2_err = out[6],
      gencorr = out[3],
      gencorr_err = out[4],
      elapsed_hours = out[7]
    )
    print(paste("Read genCorr results for",pheno,"x",CRF))
  }
}



### Appends fields name to REML and genCorr tables
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
envLM_PXS_tbl <- envLM_PXS_tbl %>% left_join(ukb_dict, by="field")
genCorr_CRF_tbl <- genCorr_CRF_tbl %>% left_join(ukb_dict, by=c("pheno_field"="field")) %>% rename(pheno_fieldname = fieldname)
genCorr_CRF_tbl <- genCorr_CRF_tbl %>% left_join(ukb_dict, by=c("CRF_field"="field")) %>% rename(CRF_fieldname = fieldname)

# Saves tables
loc_out <- paste0(dir_scratch, "REML_PXS_results.txt")
write.table(REML_PXS_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

loc_out <- paste0(dir_scratch, "REML_exposures_results.txt")
write.table(REML_expo_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

loc_out <- paste0(dir_scratch, "envLM_PXS_results.txt")
write.table(envLM_PXS_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

loc_out <- paste0(dir_scratch, "genCorr_CRF_results.txt")
write.table(genCorr_CRF_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

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
cols_to_keep <- c("SNP","CHR","BP","BETA","P_BOLT_LMM_INF","CHISQ_BOLT_LMM_INF")

sig_SNPs_PXS <- tibble(
  SNP = as.character(),
  CHR = as.numeric(),
  BP = as.numeric(),
  BETA = as.numeric(),
  P = as.numeric()
)
sig_SNPs_expo <-  sig_SNPs_PXS

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
  
  sig_SNPs_PXS <- sig_SNPs_PXS %>% add_row(
    LMM %>% select(SNP,CHR,BP,BETA,P) %>% filter(P < 5E-8)
  )
  
  N <- nrow(LMM)
  bonferroni <- 0.05 / N
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
  
  sig_SNPs_expo <- sig_SNPs_expo %>% add_row(
    LMM %>% select(SNP,CHR,BP,BETA,P) %>% filter(P < 5E-8)
  )
  
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

loc_out <- paste0(dir_scratch, "LMM_PXS_sig_SNPs.txt")
write.table(sig_SNPs_PXS, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

loc_out <- paste0(dir_scratch, "LMM_expo_sig_SNPs.txt")
write.table(sig_SNPs_expo, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

loc_out <- paste0(dir_scratch, "LMM_PXS_results.txt")
write.table(LMM_PXS_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

loc_out <- paste0(dir_scratch, "LMM_exposures_results.txt")
write.table(LMM_expo_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)
