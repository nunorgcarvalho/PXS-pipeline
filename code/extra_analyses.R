library(tidyverse)
library(data.table)
source('code/paths.R')

dir_FUMA <- paste0(dir_scratch,"FUMA_results/")
# Shared code ####
# used for plotting purposes
shortnames <- as_tibble(fread(paste0(dir_script,"../input_data/exposures_shortnames.csv")))
# loads job ID table to match jobs to fields
fields <- as_tibble(fread(paste0(dir_scratch,"fields_tbl.txt"))) %>%
  add_row(term="PXS_T2D", traitname="PXS for Type 2 Diabetes onset") %>%
  left_join(shortnames, by="term")
jobID_tbl <- as_tibble(fread(paste0(dir_FUMA,"jobid_name.csv"))) %>%
  mutate(term = substring(JobName,5))
jobID_tbl[jobID_tbl$JobName=="PXS NC", "term"] <- "PXS_T2D"
jobID_tbl <- jobID_tbl %>%
  left_join(fields %>% select(term, traitname, shortname), by="term") %>%
  select(-V4)


# GTEx BMIadj vs not ####
GTEx_general_BMIadj <- as_tibble(fread(paste0(dir_FUMA,"FUMA_job265590/magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out"), skip="VARIABLE"))
GTEx_general_NOadj <- as_tibble(fread(paste0(dir_FUMA,"FUMA_job228821/magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out"), skip="VARIABLE"))

GTEx_general <- GTEx_general_NOadj %>% select(VARIABLE, P_NOadj = P) %>%
  left_join(GTEx_general_BMIadj %>% select(VARIABLE, P_BMIadj = P), by="VARIABLE")
ggplot(GTEx_general, aes(x=-log10(P_NOadj), y=-log10(P_BMIadj))) +
  geom_point() +
  geom_abline()

# FTO SNPs ####

## extracts FUMA SNP results ####
for (i in 1:(nrow(jobID_tbl))) {
  jobID <- jobID_tbl$JobID[i]
  term <- jobID_tbl$term[i]
  shortname <- jobID_tbl$shortname[i]
  print(paste(i,jobID,term,shortname))
  
  loc_file <- paste0(dir_FUMA,"FUMA_job",jobID,"/snps.txt")
  job_data <- as_tibble(fread(loc_file, skip="uniqID"))
  
  if (i==1) {
    SNPs <- job_data %>% select(uniqID, rsID, chr, pos, gwasP, IndSigSNP, nearestGene) %>%
      mutate(term = term, shortname = shortname)
  } else {
    SNPs <- SNPs %>%
      add_row(job_data %>% select(uniqID, rsID, chr, pos, gwasP, IndSigSNP, nearestGene) %>%
      mutate(term = term, shortname = shortname) )
  }
}
FTO_terms <- (SNPs %>% filter(nearestGene == "FTO") )$term %>% unique()
FTO_SNPs <- (SNPs %>% filter(nearestGene == "FTO") )
FTO_rsIDs <- FTO_SNPs$rsID %>% unique()

window <- 50000
LMM_PXS_T2D <- as_tibble(fread("scratch/T2D/LMM_PXS_T2D_bgen.txt"))
LMM_PXS_T2D_FTO <- LMM_PXS_T2D %>%
  filter(BP >= min(FTO_SNPs$pos) - window , BP <= max(FTO_SNPs$pos) + window,
         CHR == median(FTO_SNPs$chr)) %>%
  mutate(FTO_SNP = (BP >= min(FTO_SNPs$pos) & BP <= max(FTO_SNPs$pos)) )
  #mutate(FTO_SNP = SNP %in% FTO_rsIDs )
mean_log10P <- LMM_PXS_T2D_FTO %>% group_by(FTO_SNP) %>%
  summarize(mean_log10P = mean(-log10(P_BOLT_LMM_INF)))
tt1 <- t.test(-log10(LMM_PXS_T2D_FTO[LMM_PXS_T2D_FTO$FTO_SNP,"P_BOLT_LMM_INF"]) ,
       -log10(LMM_PXS_T2D_FTO[!LMM_PXS_T2D_FTO$FTO_SNP,"P_BOLT_LMM_INF"]) )

mean_log10P <- tibble(FTO_SNP = c(F,T,F),
                      mean_log10P = c(mean(-log10( (LMM_PXS_T2D_FTO%>%filter(BP < min(FTO_SNPs$pos)))$P_BOLT_LMM_INF)),
                                      mean(-log10( (LMM_PXS_T2D_FTO%>%filter( between(BP, min(FTO_SNPs$pos), max(FTO_SNPs$pos)) ))$P_BOLT_LMM_INF)),
                                      mean(-log10( (LMM_PXS_T2D_FTO%>%filter(BP > max(FTO_SNPs$pos)))$P_BOLT_LMM_INF))
                                      ),
                      BP_min = c(min(FTO_SNPs$pos)-window, min(FTO_SNPs$pos), max(FTO_SNPs$pos)+1),
                      BP_max = c(min(FTO_SNPs$pos)-1, max(FTO_SNPs$pos), max(FTO_SNPs$pos)+window)
                      )

## mini Manhattan plot ####
gg <- ggplot(LMM_PXS_T2D_FTO, aes(x=BP, y = -log10(P_BOLT_LMM_INF))) +
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
  labs(title = "FTO SNPs don't exceed genome-wide threshold in PXS-T2D but have enriched association",
       subtitle = paste0("50kb region shown on either side of FTO SNPs. Mean -log10(P) bars shown. ",
                         "t = ", round(tt1$statistic,1), ". p = ", formatC(tt1$p.value,digits=2, format="E")) ) +
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
loc_out <- paste0(dir_Manhattan, "Zoomed_Manhattan_PXS_T2D")
ggsave(paste0(loc_out,".png"), gg, width=180, height = 120, units="mm", dpi=900)

# gvc ~ -log10(P) ####
LMM_PXS_T2D_FTO <- LMM_PXS_T2D_FTO %>%
  mutate( gvc = 2 * (BETA)^2 * (A1FREQ) * (1-A1FREQ) )
ggplot(LMM_PXS_T2D_FTO, aes(x = gvc, y = -log10(P_BOLT_LMM_INF))) +
  geom_point() +
  geom_abline(slope=1)


# gvc check ####
library(rjson)
AMP_T2D_genes <- as_tibble(fread("input_data/AMP_T2D_gene_table.csv"))
AMP_T2D_genes_sig <- AMP_T2D_genes %>% filter(minP < 2.5E-6)
PXS_T2D_genes <- as_tibble(fread("scratch/FUMA_results/FUMA_job228821/genes.txt"))
PXS_T2D_genes_sig <- PXS_T2D_genes %>% filter(minGwasP < 2.5E-6)

AMP_T2D_SNPs_json <- fromJSON(file="input_data/AMP_T2D_SNP_table.json")
AMP_T2D_SNPs <- map_dfr(AMP_T2D_SNPs_json, as_tibble) %>%
  select(varId, pValue, beta, maf, nearest) %>% distinct() %>%
  mutate( gvc = 2 * (beta)^2 * (maf) * (1-maf),
          PXS_T2D_sig_gene = nearest %in% PXS_T2D_genes_sig$symbol,
          AMP_T2D_sig_gene = nearest %in% AMP_T2D_genes_sig$gene,
          both_sig_gene = PXS_T2D_sig_gene & AMP_T2D_sig_gene) %>%
  filter(AMP_T2D_sig_gene)

ggplot(AMP_T2D_SNPs, aes(x=gvc)) +
  geom_density(aes(color = both_sig_gene))
wilcox.test(AMP_T2D_SNPs[AMP_T2D_SNPs$both_sig_gene,]$gvc,
       AMP_T2D_SNPs[!AMP_T2D_SNPs$both_sig_gene,]$gvc)
#
joint_sig_genes <- AMP_T2D_genes_sig %>% select(gene, AMP_T2D_minP = minP) %>%
  left_join(PXS_T2D_genes_sig %>% select(gene=symbol, PXS_T2D_minP = minGwasP), by="gene") %>%
  mutate(PXS_T2D_sig = !is.na(PXS_T2D_minP),
         AMP_T2D_log10P = -log10(AMP_T2D_minP))
ggplot(joint_sig_genes, aes(x=AMP_T2D_log10P)) +
  geom_density(aes(color = PXS_T2D_sig))
wilcox.test(joint_sig_genes[joint_sig_genes$PXS_T2D_sig,]$AMP_T2D_log10P,
            joint_sig_genes[!joint_sig_genes$PXS_T2D_sig,]$AMP_T2D_log10P)
