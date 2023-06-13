# Helper functions used in other scripts. No need to run independently
library(tidyverse)
library(ggrepel)


# GWAS-related functions
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
  max_logP <- ceiling(max(-log10(c(manhattan_tbl$P,bonferroni))))
  
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
                           h2_text," :: Lambda = ", round(lambda,3), " :: Lambda-adjusted P-values"))
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

# ANOVA-related functions
calculate_ANOVA <- function(field, loc_out, ANOVA_tbl, AC_means, remove_negs=FALSE, plot=FALSE) {
  
  ANOVA_tibble <- pheno %>% select(assessment_center,all_of(field)) %>%
    dplyr::rename(field_value = all_of(field)) %>%
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
  
  
  if (plot) {
    gg<-ggplot(ANOVA_tibble, aes(x=assessment_center,y=field_value)) +
      geom_boxplot() +
      geom_label(data=ANOVA_tibble_summary, aes(y=med_field_value, label=round(med_field_value,2)),
                 label.padding = unit(0.1, "lines")) +
      coord_flip() +
      xlab("Assessment Center") +
      ylab(field) +
      labs(title = paste("Boxplots of",field,"per Assessment Center"),
           subtitle = paste0("ANOVA: F = ",round(f_stat,3),", p-value = ",round(p_value,3)))
    ggsave(loc_out, gg, width = 3600, height = 2700, units = "px")
  }
  
  ANOVA_tibble_summary2 <- ANOVA_tibble_summary
  colnames(ANOVA_tibble_summary2)[c(2,3)] <- paste0(c("mean_","median_"),field)
  AC_means <- AC_means %>% left_join(ANOVA_tibble_summary2[,c(1:2)], by="assessment_center")
  
  out_list <- list(ANOVA_tbl, AC_means, ANOVA_tibble, ANOVA_tibble_summary)
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
    dplyr::rename(fieldname1 = fieldname) %>%
    left_join(ukb_dict, by=c("field2"="field")) %>%
    dplyr::rename(fieldname2 = fieldname) %>%
    arrange(p_adj)
  AC_corrs
}
