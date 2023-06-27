# Libraries and paths ####
library(tidyverse)
library(data.table)
source('paths.R')
source("old_code/helper_functions.R")

exposures_list <- readLines(paste0(dir_script,"../input_data/exposures.txt"))
shortnames <- as_tibble(fread(paste0(dir_script,"../input_data/exposures_shortnames.csv")))

# shared code ####
LMM_outs <- tibble(
    path = paste0(dir_scratch,"exposures/",exposures_list,"/LMM_",exposures_list,"_bgen.txt"),
    term = exposures_list
  ) %>% add_row(
    path = paste0(dir_scratch,"T2D/LMM_PXS_T2D_bgen.txt"),
    term = "PXS_T2D"
  ) %>%
  left_join(shortnames, by="term")

# Manhattan Plot settings
cols_to_keep <- c("SNP","CHR","BP","BETA","P_BOLT_LMM_INF","CHISQ_BOLT_LMM_INF")
colors <- c("#5BC6DD","#10404A")
minP <- 1E-30 # minimum bound for P

dir_Manhattan <- paste0(dir_script,"../final_results/Manhattan_plots/")
dir.create(dir_Manhattan, showWarnings = FALSE)

# function that filters down a large GWAS summary file to be manageable to plot
downscale_sf <- function(sf, col_P = "P_BOLT_LMM_INF", base = 4) {
  for (i in c(-1,-2,-3)) {
    top <- 10**(i+1)
    bottom <- 10**(i)
    downscale <- base**(4+i)
    allowed_i <- seq(1,nrow(sf),downscale)
    sf <- sf[c(which(!between(sf[[col_P]], bottom, top)), allowed_i),]
  }
  return(sf)
}

GWAS_summary_tbl <- tibble(term = as.character(),
                           n_SNPs = as.numeric(),
                           n_sig_SNPs = as.numeric(),
                           lambda_BOLT = as.numeric()
                           )

# Loops through and prints Manhattan plot ####
for (i in 1:nrow(LMM_outs)) {
  path <- LMM_outs$path[i]
  shortname <- LMM_outs$shortname[i]
  term <- LMM_outs$term[i]
  
  print(paste0(shortname, " :: Loading full summary file"))
  sf <- as_tibble(fread(path, select = cols_to_keep))
  # determines the cumulative base pair position of each chromosome
  manhattan_tbl <- sf %>%
    group_by(CHR) %>%
    summarize(chr_length = max(BP)) %>%
    mutate(total_bp = cumsum(as.numeric(chr_length)) - chr_length) %>%
    select(-chr_length)
  
  N <- nrow(sf) # 19400443
  #bonferroni <- 0.05 / N # very conservative bonferroni
  bonferroni <- 5E-8
  
  # gets lambda
  lambda <- median(sf$CHISQ_BOLT_LMM_INF) / qchisq(0.5,1)
  
  print(paste0(shortname, " :: Downscaling summary file"))
  # reduces the size of the plotted points to 
  sf <- sf %>% downscale_sf(base=5)
  
  # determines the cumulative base pair position of SNPs
  manhattan_tbl <- manhattan_tbl %>%
    left_join(sf, ., by="CHR") %>%
    arrange(CHR, BP) %>%
    mutate(BP_cum = total_bp + BP)
  # gets central base pair position within each chromosome
  axes_tbl <- manhattan_tbl %>%
    group_by(CHR) %>%
    summarize( center_bp = mean(c(max(BP_cum),min(BP_cum))) )
  # gets list of top SNPs
  sig_SNPs <- manhattan_tbl %>% filter(P_BOLT_LMM_INF < bonferroni) %>%
    group_by(CHR) %>% filter(P_BOLT_LMM_INF == min(P_BOLT_LMM_INF)) %>%
    unique()
  
  # sets hard limit for P values plotted
  manhattan_tbl[manhattan_tbl$P_BOLT_LMM_INF < minP,"P_BOLT_LMM_INF"] <- minP
  sig_SNPs[sig_SNPs$P_BOLT_LMM_INF < minP,"P_BOLT_LMM_INF"] <- minP
  max_logP <- ceiling(max(-log10(c(manhattan_tbl$P_BOLT_LMM_INF,bonferroni))))
  n_sig_SNPs <- nrow(manhattan_tbl %>% filter(P_BOLT_LMM_INF<bonferroni))
  
  subtitle <- paste0("Number of SNPs = ", round(N / 1E6,1), "M :: ",
                     "Number of significant SNPs = ", n_sig_SNPs,
                     " :: Lambda = ", round(lambda,2),
                     " :: Lead SNP per chromosome labeled")
  if (max_logP == -log10(minP)) {subtitle <- paste0(subtitle," :: -log10(P) capped at 30")}
  
  print(paste0(shortname, " :: Printing Manhattan plot"))
  gg<-ggplot(manhattan_tbl, aes(x = BP_cum, y = -log10(P_BOLT_LMM_INF))) +
    geom_abline(slope=0,intercept = -log10(bonferroni), color = "red") +
    geom_point(aes(color=as.factor(CHR)), alpha = 0.8, size = 0.6) +
    # geom_text_repel(data=sig_SNPs, aes(label=SNP,x=BP_cum,y=-log10(P_BOLT_LMM_INF)),
    #                 nudge_y = 0.15, seed=1, size = 0.75) +
    geom_text(data=sig_SNPs, aes(label=SNP,x=BP_cum,y=-log10(P_BOLT_LMM_INF)),
                    nudge_y = 0.3, size = 1.2) +
    scale_color_manual(values = rep(colors, 22)) +
    scale_x_continuous(label = axes_tbl$CHR, breaks = axes_tbl$center_bp, expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0:max_logP), limits = c(0,max_logP+0.5)) +
    theme_light() +
    xlab("Chromosome") +
    ylab(bquote(-log[10](P))) +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size=7),
      plot.subtitle = element_text(size=5),
      axis.title = element_text(size=6),
      axis.text = element_text(size=5)
    ) +
    labs(title = paste0("Manhattan plot for '", shortname,"'"),
         subtitle = subtitle,
         )
  #gg
  
  print(paste0(shortname, " :: Saving Manhattan plot"))
  loc_out <- paste0(dir_Manhattan, "Manhattan_",term)
  ggsave(paste0(loc_out,".png"), gg, width=180, height = 120, units="mm", dpi=900)
  if (term == "PXS_T2D") {
    ggsave(paste0(loc_out,".pdf"), gg, width=180, height = 120, units="mm", dpi=900)
  }
  
  # adds to summary table
  GWAS_summary_tbl <- GWAS_summary_tbl %>% add_row(
    term = term, n_SNPs = N, n_sig_SNPs = n_sig_SNPs, lambda_BOLT = lambda
  )
  print(paste(term, N, n_sig_SNPs, lambda))
}

# saves summary table
loc_out <- paste0(dir_results,"GWAS_summary_tbl.txt")
fwrite(GWAS_summary_tbl, loc_out, sep="\t")
