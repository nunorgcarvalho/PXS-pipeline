# Helper functions used in other scripts. No need to run independently
library(tidyverse)
library(ggrepel)


# functions
plot_Manhattan <- function(LMM, fieldname, N=-1, lambda=-1, h2=-1) {
  ### MY OWN MANHATTAN FUNCTION WITH HELP FROM https://r-graph-gallery.com/101_Manhattan_plot.html ###
  
  if (N==-1) {N <- nrow(LMM)}
  bonferroni <- 0.05 / N
  if (lambda==-1) {lambda <- get_lambda(LMM$CHISQ_BOLT_LMM_INF)}
  if (h2==-1) {
    h2_text <- ""
  } else {
    h2_text <- paste0(" :: h2 = ",round(h2,3))
  }
  
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
  
  colors <- c("#5BC6DD","#10404A")
  gg<-ggplot(manhattan_tbl, aes(x = BP_cum, y = -log10(P))) +
    geom_abline(slope=0,intercept = -log10(bonferroni), color = "red") +
    geom_point(aes(color=as.factor(CHR)), alpha = 0.8, size = 1.3) +
    geom_text_repel(data=sig_SNPs, aes(label=SNP,x=BP_cum,y=-log10(P)), nudge_y = 0.15, seed=1) +
    scale_color_manual(values = rep(colors, 22)) +
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
                           h2_text," :: Lambda = ", round(lambda,3)))
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