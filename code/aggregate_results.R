# Need at least 8GB of RAM to load full LMM data later
# Libraries and directories ####
library(tidyverse)
library(data.table)
# Functions #
source("paths.R")

# makes directory for figures
dir_figs <- paste0(dir_scratch,"figures/")
dir.create(dir_figs, showWarnings = FALSE)
dir_results <- paste0(dir_script,"../final_results/")

# paths
loc_phenolist <- paste0(dir_script,"../input_data/phenotypes.txt")
pheno_list <- readLines(loc_phenolist)
loc_expolist <- paste0(dir_script,"../input_data/exposures.txt")
exposures_list <- readLines(loc_expolist)
MAGIC_traits <- c("2hGlu","FG","FI","HbA1c")
# used for plotting purposes
shortnames <- as_tibble(fread("../input_data/exposures_shortnames.csv"))
# field information
fields <- as_tibble(fread(paste0(dir_scratch,"fields_tbl.txt"))) %>%
  left_join(shortnames, by="term") %>%
  add_row(term = c("T2D","T2D_onset","PXS_T2D", "T2D_all"), use_type="phenotype",
          shortname = c("T2D onset","T2D onset","PXS for T2D onset","Diabetes (all)")) %>%
  add_row(term = MAGIC_traits, use_type="MAGIC",
          shortname = c("2-hour Glucose","Fasting Glucose","Fasting Insulin","Glycated haemoglobin (HbA1c)"))
#manually shortens long-named traits
fields$shortname[fields$use_type=="CRF"] <- fields$traitname[fields$use_type=="CRF"]
fields[fields$term=="f4080","shortname"] <- "Systolic blood pressure"

# BOLT-REML results ####
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

## compiles list of gencorr.out results to read ####
# includes both heritability and genCorr
gencorr_logs <- tibble(path=as.character(),
                       term1=as.character(), term2=as.character())
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
## actually reads each listed result ####
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
  left_join(fields %>% select(term, shortname), by=c("pheno1_term"="term")) %>% dplyr::rename(pheno1_shortname = shortname) %>%
  left_join(fields %>% select(term, shortname), by=c("pheno2_term"="term")) %>% dplyr::rename(pheno2_shortname = shortname)

# Compares CRF gencorr with PXS_T2D vs T2D
genCorr_REML_tbl %>%
  filter(pheno1_term %in% c("PXS_T2D","T2D")) %>%
  select(pheno1_term, pheno2_shortname, gencorr, gencorr_err) %>%
  pivot_wider(names_from = pheno1_term,
              values_from = c(gencorr, gencorr_err),
              id_cols = pheno2_shortname) %>%
  arrange(-abs(gencorr_PXS_T2D))

# Compares CRF gencorr with T2D vs T2D_all
genCorr_REML_tbl %>%
  filter(pheno1_term %in% c("T2D","T2D_all")) %>%
  select(pheno1_term, pheno2_shortname, gencorr, gencorr_err) %>%
  pivot_wider(names_from = pheno1_term,
              values_from = c(gencorr, gencorr_err),
              id_cols = pheno2_shortname) %>%
  arrange(-abs(gencorr_T2D))

# Saves tables
loc_out <- paste0(dir_scratch, "genCorr_REML_results.txt")
fwrite(genCorr_REML_tbl, loc_out, sep="\t")

# prints out mean absolute gencorr for each CRF
genCorr_REML_tbl %>% filter(pheno1_term %in% exposures_list) %>%
  group_by(pheno2_term,pheno2_shortname) %>%
  summarize(h2 = mean(h2g2, na.rm=TRUE),
            mean_abs_gencorr= mean(abs(gencorr), na.rm=TRUE)) %>%
  arrange(-mean_abs_gencorr)
# prints out mean absolute gencorr for each exposure
genCorr_REML_tbl %>% filter(pheno1_term %in% exposures_list) %>%
  group_by(pheno1_term,pheno1_shortname) %>%
  summarize(h2 = mean(h2g1, na.rm=TRUE),
            mean_abs_gencorr= mean(abs(gencorr), na.rm=TRUE)) %>%
  arrange(-mean_abs_gencorr) %>% print(n=Inf)
# prints out gencorrs for just PXS_T2D x CRFs, and mean |rg|
genCorr_REML_tbl%>%filter(pheno1_term=="PXS_T2D") %>% arrange(-abs(gencorr))
(genCorr_REML_tbl%>%filter(pheno1_term=="PXS_T2D"))$gencorr %>% abs() %>% mean()

## Gencorr Visualization ####

# gets order of exposures from decreasing p-value
coeffs <- as_tibble(fread(paste0(dir_script,"../input_data/PXS_coefficients.txt")))
expo_order <- factor(c(shortnames$shortname[1],
                       (coeffs %>% filter(term %in% exposures_list) %>% 
                          left_join(fields%>%select(term,shortname), by="term") %>%
                          arrange(p.value))$shortname) )
CRF_order <- factor((genCorr_REML_tbl%>%filter(pheno1_term=="PXS_T2D") %>% arrange(-abs(gencorr)))$pheno2_shortname)
# calculates Bonferroni-corrected confidence intervals
b1 <- qnorm(1 - (0.025 / nrow(genCorr_REML_tbl %>% filter(pheno1_term %in% exposures_list))))
# Visualizes exposure x CRF gencorrs
ggplot(genCorr_REML_tbl %>% filter(pheno1_shortname %in% expo_order),
       aes(x=factor(pheno2_shortname,levels=CRF_order),
           y=fct_rev(factor(pheno1_shortname,levels=expo_order)),
           fill=gencorr)) +
  geom_tile() +
  geom_text(data=genCorr_REML_tbl %>%
              filter(pheno1_shortname %in% expo_order, abs(gencorr)-b1*gencorr_err > 0),
            aes(label = formatC(signif(gencorr,digits=2), digits=2,format="fg", flag="#"))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = length(expo_order)-0.5, ymax = length(expo_order)+0.5),
            color = "black", fill = NA) +
  scale_fill_gradient2(low="red",mid="white", high="green", midpoint = 0) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=20), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=20), expand = c(0, 0)) +
  xlab("Clinical Risk Factor (CRF)") +
  ylab("Behavior") +
  labs(title = "Genetic Correlations between behaviors + PXS-T2D and CRFs",
       subtitle = "Only significant genetic correlations shown (p < 0.05 / 150).\n
       Row for PXS-T2D highlighted. Behaviors in descending order of T2D onset association",
       fill = "Genetic Correlation") +
  theme_light()
# saves to system
loc_fig <- paste0(dir_figs,"genCorrs_PXS-expos_CRFs.png")
ggsave(loc_fig,width=4000, height=3000, units="px")

# looks at heritabilities of traits
genCorr_REML_tbl %>% group_by(pheno1_term,pheno1_shortname) %>%
  filter(h2g1_err == min(h2g1_err)) %>%
  summarize(h2 = mean(h2g1), h2_err = mean(h2g1_err)) %>%
  mutate(h2_low = h2 - 1.96*h2_err, h2_upp = h2 + 1.96*h2_err) %>%
  arrange(-h2) %>% print(n=Inf)

# LDsc results ####
## genCorr results from LDsc ####
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
### makes list of .log files to read ####
terms1 <- c("PXS_T2D","PXS_T2D","T2D_onset")
terms2 <- c("T2D_onset","T2D_all","T2D_all")

ldsc_logs <- tibble(
  path=paste0(dir_scratch,"MAGIC/MAGIC_",MAGIC_traits,"_PXS_rg.log"),
  term1="PXS_T2D", term2=MAGIC_traits
) %>% add_row(
  path=paste0(dir_scratch,"T2D/LDsc_gencorr_",terms1,"_",terms2,".log"),
  term1=terms1, term2=terms2
)
for (MAGIC_trait in MAGIC_traits) {
  ldsc_logs <- ldsc_logs %>% add_row(
    path = paste0(dir_scratch,"MAGIC/MAGIC_",MAGIC_trait,"_",exposures_list,"_rg.log"),
    term1=exposures_list, term2=MAGIC_trait
  )
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
### loops through each MAGIC trait ####
for (i in 1:nrow(ldsc_logs)) {
  # extracts ldsc gencorr data
  loc_gencorr_log <- ldsc_logs$path[i]
  if (!file.exists(loc_gencorr_log)) {next}
  gencorr_log <- readLines(loc_gencorr_log)
  out <- extract_LDsc_genCorr(gencorr_log)
  # adds to table
  genCorr_LDsc_tbl <- genCorr_LDsc_tbl %>% add_row(
    pheno1_term = ldsc_logs$term1[i], pheno2_term = ldsc_logs$term2[i],
    pheno1_h2=out[1], pheno1_h2_se=out[2], pheno1_lambda=out[3],
    pheno2_h2=out[4], pheno2_h2_se=out[5], pheno2_lambda=out[6],
    gencorr=out[7], gencorr_se=out[8], gencorr_Z=out[9], gencorr_P=out[10]
  )
}
# declares columns and uses standard z=1.96 for 95% confidence
genCorr_LDsc_tbl$gencorr_P_adj <- genCorr_LDsc_tbl$gencorr_P
genCorr_LDsc_tbl$gencorr_low <- genCorr_LDsc_tbl$gencorr - 1.96 * genCorr_LDsc_tbl$gencorr_se
genCorr_LDsc_tbl$gencorr_upp <- genCorr_LDsc_tbl$gencorr + 1.96 * genCorr_LDsc_tbl$gencorr_se
# p-adjusts gencorr(PXS_T2D, MAGIC)
# coefficient used for bonferroni-corrected confidence intervals
i_PXS_MAGIC <- which(genCorr_LDsc_tbl$pheno2_term %in% MAGIC_traits & genCorr_LDsc_tbl$pheno1_term =="PXS_T2D")
genCorr_LDsc_tbl$gencorr_P_adj[i_PXS_MAGIC] <- p.adjust(genCorr_LDsc_tbl$gencorr_P[i_PXS_MAGIC],method="bonferroni")
b2 <- qnorm(1 - (0.025 / length(i_PXS_MAGIC)))
genCorr_LDsc_tbl$gencorr_low[i_PXS_MAGIC] <- genCorr_LDsc_tbl$gencorr[i_PXS_MAGIC] - b2 * genCorr_LDsc_tbl$gencorr_se[i_PXS_MAGIC]
genCorr_LDsc_tbl$gencorr_upp[i_PXS_MAGIC] <- genCorr_LDsc_tbl$gencorr[i_PXS_MAGIC] + b2 * genCorr_LDsc_tbl$gencorr_se[i_PXS_MAGIC]
# p-adjusts gencorr(exposures, MAGIC)
i_expos_MAGIC <- which(genCorr_LDsc_tbl$pheno2_term %in% MAGIC_traits & genCorr_LDsc_tbl$pheno1_term %in% exposures_list)
genCorr_LDsc_tbl$gencorr_P_adj[i_expos_MAGIC] <- p.adjust(genCorr_LDsc_tbl$gencorr_P[i_expos_MAGIC],method="bonferroni")
b3 <- qnorm(1 - (0.025 / length(i_expos_MAGIC)))
genCorr_LDsc_tbl$gencorr_low[i_expos_MAGIC] <- genCorr_LDsc_tbl$gencorr[i_expos_MAGIC] - b2 * genCorr_LDsc_tbl$gencorr_se[i_expos_MAGIC]
genCorr_LDsc_tbl$gencorr_upp[i_expos_MAGIC] <- genCorr_LDsc_tbl$gencorr[i_expos_MAGIC] + b2 * genCorr_LDsc_tbl$gencorr_se[i_expos_MAGIC]

# determines significance
genCorr_LDsc_tbl$significant <- genCorr_LDsc_tbl$gencorr_P_adj < 0.05

# Adds named terms to table
genCorr_LDsc_tbl <- genCorr_LDsc_tbl %>%
  left_join(fields %>% select(term, shortname), by=c("pheno1_term"="term")) %>% dplyr::rename(pheno1_shortname = shortname) %>%
  left_join(fields %>% select(term, shortname), by=c("pheno2_term"="term")) %>% dplyr::rename(pheno2_shortname = shortname)

# Views significant gencorrs
genCorr_LDsc_tbl %>% filter(significant) %>%
  select(-contains("h2"), -gencorr_Z, -gencorr_P, - ends_with("lambda"), -ends_with("term"))

### Visualizes exposure x CRF gencorrs ####
ggplot(genCorr_LDsc_tbl[c(i_PXS_MAGIC,i_expos_MAGIC),],
       aes(x=pheno2_shortname,#factor(pheno2_shortname,levels=CRF_order),
           y=fct_rev(factor(pheno1_shortname,levels=expo_order)),
           fill=gencorr)) +
  geom_tile() +
  geom_text(data=genCorr_LDsc_tbl[c(i_PXS_MAGIC,i_expos_MAGIC),] %>% filter(significant),
            aes(label = formatC(signif(gencorr,digits=2), digits=2,format="fg", flag="#"))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = length(expo_order)-0.5, ymax = length(expo_order)+0.5),
            color = "black", fill = NA) +
  scale_fill_gradient2(low="red",mid="white", high="green", midpoint = 0) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=20), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=20), expand = c(0, 0)) +
  xlab("MAGIC Phenotype") +
  ylab("Behavior") +
  labs(title = "Genetic Correlations between behaviors + PXS-T2D and MAGIC phenotypes",
       subtitle = "Only significant genetic correlations shown (p < 0.05 / 100).\nRow for PXS-T2D highlighted. Behaviors in descending order of T2D onset association",
       fill = "Genetic Correlation") +
  theme_light()
# saves to system
loc_fig <- paste0(dir_figs,"genCorrs_PXS-expos_MAGIC.png")
ggsave(loc_fig,width=4000, height=3000, units="px")











## h2 results from LDsc ####
# function for extracting data from ldsc h2 log file
extract_LDsc_h2 <- function(h2_log) {
  N <- length(h2_log)
  
  # extract h2
  line <- h2_log[N-6]
  line_split <- str_split(line, ": ")[[1]][2]
  h2 <- as.numeric(str_split(line_split, " ")[[1]][1])
  h2_se <- parse_number(str_split(line_split, " ")[[1]][2])
  # extract lambda
  line <- h2_log[N-5]
  line_split <- str_split(line, ": ")[[1]][2]
  lambda <- as.numeric(line_split)
  # extract mean_chi2
  line <- h2_log[N-4]
  line_split <- str_split(line, ": ")[[1]][2]
  mean_chi2 <- as.numeric(line_split)
  # extract intercept
  line <- h2_log[N-3]
  line_split <- str_split(line, ": ")[[1]][2]
  intercept <- as.numeric(str_split(line_split, " ")[[1]][1])
  intercept_se <- parse_number(str_split(line_split, " ")[[1]][2])
  # extract ratio
  line <- h2_log[N-2]
  line_split <- str_split(line, ": ")[[1]][2]
  ratio <- as.numeric(str_split(line_split, " ")[[1]][1])
  ratio_se <- parse_number(str_split(line_split, " ")[[1]][2])
  
  out <- c(h2, h2_se, lambda, mean_chi2, intercept, intercept_se, ratio, ratio_se)
  out
}
### makes list of ldsc log files to go through
ldsc_h2_logs <- tibble(
    path = paste0(dir_scratch,"exposures/",exposures_list,"/ldsc_h2_",exposures_list,".log"),
    term = exposures_list
  ) %>% add_row( path = paste0(dir_scratch,"T2D/ldsc_h2_PXS_T2D.log"), term = "PXS_T2D" )

h2_LDsc_tbl <- tibble(
  term = as.character(),
  h2 = as.numeric(),
  h2_se = as.numeric(),
  lambda = as.numeric(),
  mean_chi2 = as.numeric(),
  intercept = as.numeric(),
  intercept_se = as.numeric(),
  ratio = as.numeric(),
  ratio_se = as.numeric()
)
### loops through each exposure ####
for (i in 1:nrow(ldsc_h2_logs)) {
  loc_h2_log <- ldsc_h2_logs$path[i]
  term <- ldsc_h2_logs$term[i]
  #loc_h2_log <- paste0(dir_scratch,"exposures/",expo,"/ldsc_h2_",expo,".log")
  # skips if log file is missing
  if (!file.exists(loc_h2_log)) {print(paste("Missing log file for",term)); next}
  print(paste("Reading log file for", term))
  h2_log <- readLines(loc_h2_log)
  if (length(h2_log)==0) {print("Log file is empty!"); next}
  out <- extract_LDsc_h2(h2_log)
  # adds data to table
  h2_LDsc_tbl <- h2_LDsc_tbl %>% add_row(
    term=term, h2=out[1], h2_se=out[2], lambda=out[3], mean_chi2=out[4],
    intercept=out[5], intercept_se=out[6], ratio=out[7], ratio_se=out[8]
  )
}
# appends trait names
h2_LDsc_tbl <- h2_LDsc_tbl %>% left_join(fields[,c("term","shortname")], by="term")

# compares BOLT-REML vs LDsc in h2 computation
h2_comparison <- h2_LDsc_tbl %>%
  select(term, shortname, ldsc_h2=h2,ldsc_h2_se=h2_se,ldsc_lambda=lambda, ldsc_ratio=ratio) %>%
  left_join(genCorr_REML_tbl %>% group_by(pheno1_term) %>%
              summarize(REML_h2 = mean(h2g1), REML_h2_err = mean(h2g1_err)),
            by = c("term"="pheno1_term"))
### ldsc vs REML comparison ####
# shows h2 comparisons after LDsc's GC adjustment
(lm1 <- lm(ldsc_h2 ~ REML_h2, data = h2_comparison))
ggplot(h2_comparison, aes(x=REML_h2, y=ldsc_h2)) +
  geom_abline(slope=1) +
  geom_smooth(method='lm') +
  geom_errorbar(aes(ymin = ldsc_h2-2*ldsc_h2_se, ymax = ldsc_h2+2*ldsc_h2_se)) +
  geom_errorbarh(aes(xmin = REML_h2-2*REML_h2_err, xmax = REML_h2+2*REML_h2_err)) +
  geom_point() +
  xlim(0,0.2) + ylim(0,0.2) +
  xlab("h2 estimated by BOLT-REML") +
  ylab("h2 estimated by LDsc") +
  labs(title = "Comparison of h2 estimate of 25 exposures by LDsc vs BOLT-REML",
       subtitle = paste0("LDsc_h2 ~ BOLT-REML_h2 slope = ", round(lm1$coefficients[2],3)))

# prints table
loc_out <- paste0(dir_results, "h2_ldsc_REML.txt")
fwrite(h2_comparison, loc_out, sep="\t")
