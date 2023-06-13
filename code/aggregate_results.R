# Need at least 8GB of RAM to load full LMM data later
## Libraries and directories ##
library(tidyverse)
library(data.table)
# Functions #
source("paths.R")

# makes directory for figures
dir_figs <- paste0(dir_scratch,"figures/")
dir.create(dir_figs, showWarnings = FALSE)

# paths
loc_phenolist <- paste0(dir_script,"../input_data/phenotypes.txt")
pheno_list <- readLines(loc_phenolist)
loc_expolist <- paste0(dir_script,"../input_data/exposures.txt")
exposures_list <- readLines(loc_expolist)
fields <- as_tibble(fread(paste0(dir_scratch,"fields_tbl.txt"))) %>%
  add_row(term = c("T2D","PXS_T2D", "T2D_all"),
          traitname = c("Type 2 Diabetes onset","PXS for Type 2 Diabetes onset","Type 2 Diabetes (all)"))
#manually shortens long-named traits
fields[fields$term=="f4080","traitname"] <- "Systolic blood pressure"

## Functions ##
extract_REML_genCorr <- function(genCorr) {
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

## Code ##

### REML results for genCorr PXS with CRFs
# includes both heritability and genCorr
gencorr_logs <- tibble(path=as.character(),
                           term1=as.character(), term2=as.character())
# compiles list of gencorr.out results to read
for (pheno in pheno_list) {
  loc_CRFs <- paste0(dir_scratch,pheno,"/",pheno,"_CRFs.txt")
  if (!file.exists(loc_CRFs)) {next}
  CRFs_list <- readLines(loc_CRFs)
  for (CRF in CRFs_list) {
    # reads exposures x CRFs gencorrs
    for (expo in exposures_list) {
      loc_genCorr <- paste0(dir_scratch,"exposures/",expo,"/",CRF,"_",expo,"_genCorr.out")
      if (!file.exists(loc_genCorr)) {
        print(paste("Did not find genCorr for",expo,"x",CRF))
        next
      }
      gencorr_logs <- gencorr_logs %>%
        add_row(path=loc_genCorr,term1=expo,term2=CRF)
    }
    
    # reads PXS x CRFs gencorrs
    loc_genCorr <- paste0(dir_scratch,pheno,"/PXS_",pheno,"_",CRF,"_genCorr.out")
    gencorr_logs <- gencorr_logs %>%
      add_row(path=loc_genCorr,term1=paste0("PXS_",pheno),term2=CRF)
    
    # reads T2D x CRFs gencorrs
    loc_genCorr <- paste0(dir_scratch,pheno,"/",pheno,"_",CRF,"_genCorr.out")
    gencorr_logs <- gencorr_logs %>%
      add_row(path=loc_genCorr,term1=pheno,term2=CRF)
    
    # reads T2D_all x CRFs gencorrs
    loc_genCorr <- paste0(dir_scratch,pheno,"/",pheno,"_all_",CRF,"_genCorr.out")
    gencorr_logs <- gencorr_logs %>%
      add_row(path=loc_genCorr,term1=paste0(pheno,"_all"),term2=CRF)
  }
}
genCorr_REML_tbl <- tibble(
  pheno1_term = as.character(),
  h2g1 = as.numeric(),
  h2g1_err = as.numeric(),
  pheno2_term = as.character(),
  h2g2 = as.numeric(),
  h2g2_err = as.numeric(),
  gencorr = as.numeric(),
  gencorr_err = as.numeric(),
  elapsed_hours = as.numeric()
)
# actually reads each listed result
for (i in 1:nrow(gencorr_logs)) {
  gencorr_log <- readLines(gencorr_logs$path[i])
  out <- extract_REML_genCorr(gencorr_log)
  genCorr_REML_tbl <- genCorr_REML_tbl %>% add_row(
    pheno1_term = gencorr_logs$term1[i],
    h2g1 = out[1],
    h2g1_err = out[2],
    pheno2_term = gencorr_logs$term2[i],
    h2g2 = out[5],
    h2g2_err = out[6],
    gencorr = out[3],
    gencorr_err = out[4],
    elapsed_hours = out[7]
  )
  print(paste("Read genCorr results for",gencorr_logs$term1[i],"x",gencorr_logs$term2[i]))
  
}
# Prints out missing data
print("GenCorr combinations with missing data:")
genCorr_REML_tbl %>% filter(rowSums(is.na(.)) > 0) %>% select(ends_with("term"))

# Adds named terms to table
genCorr_REML_tbl <- genCorr_REML_tbl %>%
  left_join(fields %>% select(term, traitname), by=c("pheno1_term"="term")) %>% rename(pheno1_traitname = traitname) %>%
  left_join(fields %>% select(term, traitname), by=c("pheno2_term"="term")) %>% rename(pheno2_traitname = traitname)

# Compares CRF gencorr with PXS_T2D vs T2D
genCorr_REML_tbl %>%
  filter(pheno1_term %in% c("PXS_T2D","T2D")) %>%
  select(pheno1_term, pheno2_traitname, gencorr, gencorr_err) %>%
  pivot_wider(names_from = pheno1_term,
              values_from = c(gencorr, gencorr_err),
              id_cols = pheno2_traitname) %>%
  arrange(-abs(gencorr_PXS_T2D))

# Compares CRF gencorr with T2D vs T2D_all
genCorr_REML_tbl %>%
  filter(pheno1_term %in% c("T2D","T2D_all")) %>%
  select(pheno1_term, pheno2_traitname, gencorr, gencorr_err) %>%
  pivot_wider(names_from = pheno1_term,
              values_from = c(gencorr, gencorr_err),
              id_cols = pheno2_traitname) %>%
  arrange(-abs(gencorr_T2D))

# Saves tables
loc_out <- paste0(dir_scratch, "genCorr_REML_results.txt")
write.table(genCorr_REML_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

# prints out mean absolute gencorr for each CRF
genCorr_REML_tbl %>% filter(pheno1_term %in% exposures_list) %>%
  group_by(pheno2_term,pheno2_traitname) %>%
  summarize(h2 = mean(h2g2, na.rm=TRUE),
            mean_abs_gencorr= mean(abs(gencorr), na.rm=TRUE)) %>%
  arrange(-mean_abs_gencorr)
# prints out mean absolute gencorr for each exposure
genCorr_REML_tbl %>% filter(pheno1_term %in% exposures_list) %>%
  group_by(pheno1_term,pheno1_traitname) %>%
  summarize(h2 = mean(h2g1, na.rm=TRUE),
            mean_abs_gencorr= mean(abs(gencorr), na.rm=TRUE)) %>%
  arrange(-mean_abs_gencorr) %>% print(n=Inf)
# prints out gencorrs for just PXS_T2D x CRFs, and mean |rg|
genCorr_REML_tbl%>%filter(pheno1_term=="PXS_T2D") %>% arrange(-abs(gencorr))
(genCorr_REML_tbl%>%filter(pheno1_term=="PXS_T2D"))$gencorr %>% abs() %>% mean()

# calculates Bonferroni-corrected confidence intervals
b1 <- qnorm(1 - (0.025 / nrow(genCorr_REML_tbl %>% filter(pheno1_term %in% exposures_list))))
# gets order of exposures from decreasing p-value
coeffs <- as_tibble(fread(paste0(dir_script,"../input_data/PXS_coefficients.txt")))
expo_order <- factor(c("PXS for Type 2 Diabetes onset",
                       (coeffs %>%
                          filter(term %in% exposures_list) %>% 
                          left_join(fields%>%select(term,traitname), by="term") %>%
                          arrange(p.value))$traitname)
                     )
CRF_order <- factor((genCorr_REML_tbl%>%filter(pheno1_term=="PXS_T2D") %>% arrange(-abs(gencorr)))$pheno2_traitname)
# Visualizes exposure x CRF gencorrs
ggplot(genCorr_REML_tbl %>% filter(pheno1_traitname %in% expo_order),
       aes(x=factor(pheno2_traitname,levels=CRF_order),
           y=fct_rev(factor(pheno1_traitname,levels=expo_order)),
           fill=gencorr)) +
  geom_tile() +
  geom_text(data=genCorr_REML_tbl %>%
              filter(pheno1_traitname %in% expo_order, abs(gencorr)-b1*gencorr_err > 0),
            aes(label = formatC(signif(gencorr,digits=2), digits=2,format="fg", flag="#"))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = length(expo_order)-0.5, ymax = length(expo_order)+0.5),
            color = "black", fill = NA) +
  scale_fill_gradient2(low="red",mid="white", high="green", midpoint = 0) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=20), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=40), expand = c(0, 0)) +
  xlab("Clinical Risk Factor (CRF)") +
  ylab("Behavior") +
  labs(title = "Genetic Correlations between behaviors + PXS-T2D and CRFs",
       subtitle = "Only significant genetic correlations shown (p < 0.05 / 150). Row for PXS-T2D highlighted",
       fill = "Genetic Correlation") +
  theme_light()
# saves to system
loc_fig <- paste0(dir_figs,"genCorrs_PXS-expos_CRFs.png")
ggsave(loc_fig,width=4000, height=3000, units="px")



#### genCorr results from LDsc ####
extract_LDsc_genCorr <- function(gencorr_log) {
  N <- length(gencorr_log)
  
  # extract pheno1_h2
  line <- gencorr_log[N-32]
  line_split <- str_split(line, ": ")[[1]][2]
  pheno1_h2 <- as.numeric(str_split(line_split, " ")[[1]][1])
  pheno1_h2_se <- parse_number(str_split(line_split, " ")[[1]][2])
  # extract pheno1_lambda
  line <- gencorr_log[N-31]
  line_split <- str_split(line, ": ")[[1]][2]
  pheno1_lambda <- as.numeric(line_split)
  
  
  # extract pheno2_h2
  line <- gencorr_log[N-24]
  line_split <- str_split(line, ": ")[[1]][2]
  pheno2_h2 <- as.numeric(str_split(line_split, " ")[[1]][1])
  pheno2_h2_se <- parse_number(str_split(line_split, " ")[[1]][2])
  # extract pheno1_lambda
  line <- gencorr_log[N-23]
  line_split <- str_split(line, ": ")[[1]][2]
  pheno2_lambda <- as.numeric(line_split)
  
  # extract gencorr
  line <- gencorr_log[N-10]
  line_split <- str_split(line, ": ")[[1]][2]
  gencorr <- as.numeric(str_split(line_split, " ")[[1]][1])
  gencorr_se <- parse_number(str_split(line_split, " ")[[1]][2])
  
  # extract gencorr_Z and gencorr P
  line <- gencorr_log[N-9]
  gencorr_Z <- as.numeric(str_split(line, ": ")[[1]][2])
  line <- gencorr_log[N-8]
  gencorr_P <- as.numeric(str_split(line, ": ")[[1]][2])
  
  out <- c(pheno1_h2, pheno1_h2_se, pheno1_lambda, pheno2_h2, pheno2_h2_se, pheno2_lambda,
           gencorr, gencorr_se, gencorr_Z, gencorr_P)
  out
}

genCorr_LDsc_tbl <- tibble(
  pheno1_term = as.character(),
  pheno1_h2 = as.numeric(),
  pheno1_h2_se = as.numeric(),
  pheno1_lambda = as.numeric(),
  pheno2_term = as.character(),
  pheno2_h2 = as.numeric(),
  pheno2_h2_se = as.numeric(),
  pheno2_lambda = as.numeric(),
  gencorr = as.numeric(),
  gencorr_se = as.numeric(),
  gencorr_Z = as.numeric(),
  gencorr_P = as.numeric()
)
MAGIC_traits <- c("2hGlu","FG","FI","HbA1c")
for (MAGIC_trait in MAGIC_traits) {
  loc_gencorr_log <- paste0("../scratch/MAGIC/MAGIC_",MAGIC_trait,"_PXS_rg.log")
  gencorr_log <- readLines(loc_gencorr_log)
  out <- extract_LDsc_genCorr(gencorr_log)
  
  genCorr_LDsc_tbl <- genCorr_LDsc_tbl %>% add_row(
    pheno1_term = MAGIC_trait, pheno2_term = "PXS_T2D",
    pheno1_h2=out[1], pheno1_h2_se=out[2], pheno1_lambda=out[3],
    pheno2_h2=out[4], pheno2_h2_se=out[5], pheno2_lambda=out[6],
    gencorr=out[7], gencorr_se=out[8], gencorr_Z=out[9], gencorr_P=out[10]
  )
}
b2 <- qnorm(1 - (0.025 / nrow(genCorr_LDsc_tbl)))
genCorr_LDsc_tbl$gencorr_P_adj <- p.adjust(genCorr_LDsc_tbl$gencorr_P, method="bonferroni")
genCorr_LDsc_tbl <- genCorr_LDsc_tbl %>%
  mutate(gencorr_low = gencorr - b2 * gencorr_se,
         gencorr_upp = gencorr + b2 * gencorr_se)
genCorr_LDsc_tbl %>% select(-ends_with("_se"),-ends_with("lambda"))

#### GWAS ANALYSIS ########

# # Code
# source("helper_functions.R")
# cols_to_keep <- c("SNP","CHR","BP","BETA","P_BOLT_LMM_INF","CHISQ_BOLT_LMM_INF")
# 
# sig_SNPs_PXS <- tibble(
#   SNP = as.character(),
#   CHR = as.numeric(),
#   BP = as.numeric(),
#   BETA = as.numeric(),
#   P_unadj = as.numeric(),
#   P = as.numeric(),
#   field = as.character()
# )
# sig_SNPs_expo <-  sig_SNPs_PXS
# 
# LMM_PXS_tbl <- tibble(
#   field = as.character(),
#   lambda = as.numeric(),
#   N = as.numeric(),
#   N_above_bonferroni = as.numeric(),
#   N_above_5E_8 = as.numeric()
# )
# LMM_expo_tbl <- LMM_PXS_tbl
# 
# # Loop through disease PXSs
# for (pheno in pheno_list) {
#   
#   fieldname <- paste0(pheno, " PXS")
#   loc_LMM <- paste0(dir_scratch,pheno,"/LMM_",pheno,"_bgen.txt")
#   LMM <- as_tibble(fread(loc_LMM, select=cols_to_keep)) %>%
#     rename(P_unadj = P_BOLT_LMM_INF)
#   
#   N <- nrow(LMM)
#   bonferroni <- 0.05 / N
#   lambda <- get_lambda(LMM$CHISQ_BOLT_LMM_INF)
#   
#   # adjusts for lambda inflation factor
#   LMM <- LMM %>% mutate(
#     P = exp(pchisq(CHISQ_BOLT_LMM_INF/lambda,1, lower.tail=FALSE, log.p=TRUE))
#   )
#   sig_SNPs_PXS <- sig_SNPs_PXS %>% add_row(
#     LMM %>% select(SNP,CHR,BP,BETA,P_unadj,P) %>% filter(P_unadj < 5E-8) %>% mutate(field = pheno)
#   )
#   LMM_PXS_tbl <- LMM_PXS_tbl %>%
#     add_row(
#       field = pheno,
#       lambda = lambda,
#       N = N,
#       N_above_bonferroni = nrow(LMM %>% filter(P < bonferroni)),
#       N_above_5E_8 = nrow(LMM %>% filter(P < 5E-8))
#     )
#   
#   
#   # reduces the number of points to plot significantly while trying to keep
#   # the Manhattan plot visually the same as before
#   LMM2 <- downscale_sf(LMM)
#   
#   gg <- plot_Manhattan(LMM2, fieldname, N, lambda)
#   
#   loc_out <- paste0(dir_scratch,pheno,"/LMM_",pheno,"_manhattan.png")
#   ggsave(loc_out, gg, width=3600,height=2700, units="px")
#   
#   print(paste("Saved Manhattan plot for",pheno))
# }
# # Loop through exposures
# for (expo in exposures_list) {
#   
#   fieldname <- expo
#   loc_LMM <- paste0(dir_scratch,"exposures/",expo,"/LMM_",expo,"_bgen.txt")
#   if (!file.exists(loc_LMM)) {next}
#   LMM <- as_tibble(fread(loc_LMM, select=cols_to_keep)) %>%
#     rename(P_unadj = P_BOLT_LMM_INF)
#   
#   N <- nrow(LMM)
#   bonferroni <- 0.05 / N
#   lambda <- get_lambda(LMM$CHISQ_BOLT_LMM_INF)
#   
#   # adjusts for lambda inflation factor
#   if (lambda > 1) {
#     LMM <- LMM %>% mutate(
#       P = exp(pchisq(CHISQ_BOLT_LMM_INF/lambda,1, lower.tail=FALSE, log.p=TRUE))
#     )
#   } else { LMM <- LMM %>% mutate(P = P_unadj) }
#   sig_SNPs_expo <- sig_SNPs_expo %>% add_row(
#     LMM %>% select(SNP,CHR,BP,BETA,P_unadj,P) %>% filter(P_unadj < 5E-8) %>% mutate(field = expo)
#   )
#   
#   LMM_expo_tbl <- LMM_expo_tbl %>%
#     add_row(
#       field = expo,
#       lambda = lambda,
#       N = N,
#       N_above_bonferroni = nrow(LMM %>% filter(P < bonferroni)),
#       N_above_5E_8 = nrow(LMM %>% filter(P < 5E-8))
#     )
#   
#   
#   # reduces the number of points to plot significantly while trying to keep
#   # the Manhattan plot visually the same as before
#   LMM2 <- downscale_sf(LMM)
#   
#   gg <- plot_Manhattan(LMM2, fieldname, N, lambda)
#   
#   loc_out <- paste0(dir_scratch,"exposures/",expo,"/LMM_",expo,"_manhattan.png")
#   ggsave(loc_out, gg, width=3600,height=2700, units="px")
#   
#   print(paste("Saved Manhattan plot for",expo))
# }
# 
# 
# LMM_PXS_tbl <- LMM_PXS_tbl %>% left_join(ukb_dict, by="field")
# LMM_expo_tbl <- LMM_expo_tbl %>% left_join(fields %>% select(term,fieldname, Meaning), by=c("field"="term"))
# 
# loc_out <- paste0(dir_scratch, "LMM_PXS_sig_SNPs.txt")
# write.table(sig_SNPs_PXS, loc_out, sep="\t", quote=FALSE, row.names=FALSE)
# 
# loc_out <- paste0(dir_scratch, "LMM_expo_sig_SNPs.txt")
# write.table(sig_SNPs_expo, loc_out, sep="\t", quote=FALSE, row.names=FALSE)
# 
# loc_out <- paste0(dir_scratch, "LMM_PXS_results.txt")
# write.table(LMM_PXS_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)
# 
# loc_out <- paste0(dir_scratch, "LMM_exposures_results.txt")
# write.table(LMM_expo_tbl, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

## Looking closer at envLM
# AC_tbl <- as_tibble(fread(paste0(dir_data_showcase,"Codings.tsv"),quote="")) %>%
#   filter(Coding == 10)
# envLM_tval_group <- envLM_PXS_tvals %>%
#   mutate(tvalue = abs(tvalue)) %>%
#   arrange(-tvalue)
# envLM_tval_group <- envLM_tval_group %>%
#   group_by(variable) %>%
#   summarize(tvalue = mean(tvalue),
#             n = n()) %>%
#   arrange(-tvalue)
# ACs <- c()
# for (i in 1:nrow(envLM_tval_group)) {
#   if (substr(envLM_tval_group$variable[i],1,17) == "assessment_center") {
#     code <- substr(envLM_tval_group$variable[i],18,22)
#     ACs[i] <- (AC_tbl %>% filter(Value==code))$Meaning[1]
#   } else {ACs[i] <- as.character(NA)}
# }
# envLM_tval_group$assessment_center <- ACs

