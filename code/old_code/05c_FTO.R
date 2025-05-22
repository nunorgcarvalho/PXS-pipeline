library(tidyverse)
library(data.table)
source('code/paths.R')

dir_FUMA <- paste0(dir_scratch,"FUMA_results/")
# Shared code ####

# FTO SNPs ####
list_FUMA_SNPs <- list.files(dir_FUMA, pattern='snps.txt', recursive=TRUE)
## extracts FUMA SNP results ####
for (i in 1:length(list_FUMA_SNPs)) {
  file <- list_FUMA_SNPs[i]
  loc_file <- paste0(dir_FUMA, file)
  print(paste(i,file))
  
  job_data <- as_tibble(fread(loc_file, skip="uniqID"))
  if (nrow(job_data)==0) {next}
  
  if (i==1) {
    SNPs <- job_data %>% select(uniqID, rsID, chr, pos, gwasP, IndSigSNP, nearestGene)
  } else {
    SNPs <- SNPs %>%
      add_row(job_data %>% select(uniqID, rsID, chr, pos, gwasP, IndSigSNP, nearestGene))
  }
}
FTO_SNPs <- (SNPs %>% filter(nearestGene == "FTO") )
FTO_rsIDs <- FTO_SNPs$rsID %>% unique()

window <- 250000
LMM_BRS <- as_tibble(fread(paste0(dir_scratch, "LMM/LMM.ALL.BRS-ALL-cov_bvr.bgen.txt")))
LMM_BRS_FTO <- LMM_BRS %>%
  filter(BP >= min(FTO_SNPs$pos) - window , BP <= max(FTO_SNPs$pos) + window,
         CHR == median(FTO_SNPs$chr)) %>%
  mutate(FTO_SNP = (BP >= min(FTO_SNPs$pos) & BP <= max(FTO_SNPs$pos)) )

LMM_BRS_FTO %>% arrange(P_BOLT_LMM_INF) %>% select(P_BOLT_LMM_INF, everything()) %>%
  filter(FTO_SNP) %>% View()
mean_log10P <- LMM_BRS_FTO %>% group_by(FTO_SNP) %>%
  summarize(mean_log10P = mean(-log10(P_BOLT_LMM_INF)))
tt1 <- t.test(-log10(LMM_BRS_FTO[LMM_BRS_FTO$FTO_SNP,"P_BOLT_LMM_INF"]) ,
       -log10(LMM_BRS_FTO[!LMM_BRS_FTO$FTO_SNP,"P_BOLT_LMM_INF"]) )

mean_log10P <- tibble(FTO_SNP = c(F,T,F),
                      mean_log10P = c(mean(-log10( (LMM_BRS_FTO%>%filter(BP < min(FTO_SNPs$pos)))$P_BOLT_LMM_INF)),
                                      mean(-log10( (LMM_BRS_FTO%>%filter( between(BP, min(FTO_SNPs$pos), max(FTO_SNPs$pos)) ))$P_BOLT_LMM_INF)),
                                      mean(-log10( (LMM_BRS_FTO%>%filter(BP > max(FTO_SNPs$pos)))$P_BOLT_LMM_INF))
                                      ),
                      BP_min = c(min(FTO_SNPs$pos)-window, min(FTO_SNPs$pos), max(FTO_SNPs$pos)+1),
                      BP_max = c(min(FTO_SNPs$pos)-1, max(FTO_SNPs$pos), max(FTO_SNPs$pos)+window)
                      )

LMM_BRS_FTO %>% filter(FTO_SNP, P_BOLT_LMM_INF < 1E-5)
## mini Manhattan plot ####
gg <- ggplot(LMM_BRS_FTO, aes(x=BP, y = -log10(P_BOLT_LMM_INF))) +
  geom_hline(yintercept = -log10(5E-8), color = "red") +
  geom_hline(yintercept = -log10(1E-5), color = "blue") +
  geom_point(aes(color = FTO_SNP), alpha = 0.5, size=1) +
  geom_segment(data = mean_log10P, aes(x = BP_min, xend = BP_max, color = FTO_SNP,
                                       y = mean_log10P, yend = mean_log10P),
               size=2) +
  xlab("Base pair (chromosome 16)") +
  ylab(bquote(-log[10](P))) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0:8), limits = c(0,8)) +
  scale_color_manual(values = c("#93CF1A","#557F00")) +
  # labs(title = "FTO SNPs don't exceed genome-wide threshold in PXS-T2D but have enriched association",
  #      subtitle = paste0("50kb region shown on either side of FTO SNPs. Mean -log10(P) bars shown. ",
  #                        "t = ", round(tt1$statistic,1), ". p = ", formatC(tt1$p.value,digits=2, format="E")) ) +
  theme_light() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=7),
    plot.subtitle = element_text(size=5),
    axis.title = element_text(size=6),
    axis.text = element_text(size=5)
  )
gg
dir_Manhattan <- paste0(dir_script,"../final_results/Manhattan_plots/")
loc_out <- paste0(dir_Manhattan, "Zoomed_Manhattan_BRS")
ggsave(paste0(loc_out,".png"), gg, width=180, height = 120, units="mm", dpi=900)

# gvc ~ -log10(P) ####
LMM_BRS_FTO <- LMM_BRS_FTO %>%
  mutate( gvc = 2 * (BETA)^2 * (A1FREQ) * (1-A1FREQ) )
ggplot(LMM_BRS_FTO, aes(x = gvc, y = -log10(P_BOLT_LMM_INF))) +
  geom_point() +
  geom_abline(slope=1)