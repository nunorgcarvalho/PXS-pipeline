#### Libraries and paths ####
library(tidyverse)
library(data.table)
library(ggrepel)
source('paths.R')

dir_FUMA <- paste0(dir_scratch,"FUMA_results/")
# loads job ID table to match jobs to fields
fields <- as_tibble(fread(paste0(dir_scratch,"fields_tbl.txt")))
jobID_tbl <- as_tibble(fread(paste0(dir_FUMA,"jobid_name.csv"))) %>%
  mutate(term = substring(JobName,5)) %>%
  left_join(fields %>% select(term, traitname), by="term") %>%
  select(-V4)
jobID_tbl[jobID_tbl$JobName=="PXS NC",c("term","traitname")] <- list("PXS_T2D","PXS for Type 2 Diabetes onset")

#### GTEx Data ####
# gets order of exposures from decreasing p-value
coeffs <- as_tibble(fread(paste0(dir_script,"../input_data/PXS_coefficients.txt")))
expo_order <- factor(c("PXS for Type 2 Diabetes onset",
                       (coeffs %>%
                          filter(term %in% fields[fields$use_type=="exposure",]$term) %>% 
                          left_join(fields%>%select(term,traitname), by="term") %>%
                          arrange(p.value))$traitname) )

# looks at GTEx general data
for (i in 1:(nrow(jobID_tbl))) {
  jobID <- jobID_tbl$JobID[i]
  term <- jobID_tbl$term[i]
  traitname <- jobID_tbl$traitname[i]
  print(paste(i,jobID,term,traitname))
  
  loc_file <- paste0(dir_FUMA,"FUMA_job",jobID,"/magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out")
  job_data <- as_tibble(fread(loc_file, skip="VARIABLE"))
  
  if (i==1) {GTEx_general <- job_data %>% mutate(term = term, traitname = traitname)
  } else {
    GTEx_general <- GTEx_general %>%
      add_row(job_data %>% mutate(term = term, traitname = traitname))
  }
}
# Bonferroni correction (FUMA already Bonferonni corrects within each trait)
GTEx_general$P_adjusted <- p.adjust(GTEx_general$P / length(unique(GTEx_general$VARIABLE)), method = 'bonferroni')
GTEx_general$significant_P_adj <- GTEx_general$P_adjusted < 0.05
GTEx_general %>% filter(significant_P_adj) %>% arrange(P_adjusted) %>% print(n=Inf)
# tile plot
ggplot(GTEx_general, aes(x = VARIABLE,
                         y = fct_rev(factor(traitname, levels=expo_order)))) +
  geom_tile(aes(fill=-log10(P_adjusted))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = length(expo_order)-0.5, ymax = length(expo_order)+0.5),
            color = "black", fill = NA) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=30), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=50), expand = c(0, 0)) +
  scale_fill_gradient2(low="white",mid="white", high="green", midpoint = -log10(0.05)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Tissue Category") + ylab("Trait") +
  labs(title = "Average gene expression per tissue category for each behavior",
       subtitle = "Using MAGMA with GTEx through FUMA. White: p<0.05 (Bonferroni adjusted)")

# compares brain expression to T2D_onset association
Brain_PXS_tbl <- GTEx_general %>%
  filter(VARIABLE=="Brain", term != "PXS_T2D") %>%
  left_join(coeffs %>% select(term, estimate, p.value), by="term") %>%
  rename(Brain_P = P_adjusted, PXS_P = p.value) #%>% filter(-log10(PXS_P) < 20)
(cor1 <- cor.test(rank(Brain_PXS_tbl$Brain_P), rank(Brain_PXS_tbl$PXS_P)))
ggplot(Brain_PXS_tbl, aes(x = rank(Brain_P), y = rank(PXS_P))) +
  geom_smooth(method='lm') +
  geom_point() +
  geom_text_repel(aes(label = traitname)) +
  scale_x_reverse() +
  scale_y_reverse() +
  xlab("Rank of behavior-brain tissue expression enrichment (by P-value)") +
  ylab("Rank of T2D onset association strength in PXS (by P-value)") +
  labs(title = "Relationship between T2D onset strength and brain tissue expression enrichment",
       subtitle = paste0("r = ", round(cor1$estimate,3), " :: p = ", round(cor1$p.value,3))) +
  theme_light()

GTEx_general %>% filter(term == "PXS_T2D") %>% arrange(P_adjusted) %>%
  select(VARIABLE, P_adjusted, significant_P_adj) %>% print(n=15)
GTEx_specific %>% filter(term == "PXS_T2D") %>% arrange(P_adjusted) %>%
  select(FULL_NAME, P_adjusted, significant_P_adj) %>% print(n=15)

# repeats the above correlation for all tissue types
tissue_PXS_GTEx <- GTEx_specific %>%
  filter(term != "PXS_T2D") %>%
  select(VARIABLE, P, term, traitname) %>%
  pivot_wider(
    names_from = VARIABLE,
    names_prefix = "P_",
    values_from = P
  ) %>%
  left_join(coeffs %>% select(term, estimate, p.value), by="term")
PXS_rank <- rank(tissue_PXS_GTEx$p.value)
#PXS_rank <- (tissue_PXS_GTEx$p.value)
tissue_PXS_GTEx_results <- tibble(
  tissue = as.character(), r = as.numeric(), p = as.numeric()
)
for (tissue in unique(GTEx_specific$VARIABLE)) {
  tissue_GTEx_rank <- rank(tissue_PXS_GTEx[[paste0("P_",tissue)]])
  #tissue_GTEx_rank <- (tissue_PXS_GTEx[[paste0("P_",tissue)]])
  cor_tissue <- cor.test(PXS_rank, tissue_GTEx_rank)
  print(paste(tissue, round(cor_tissue$estimate,3), round(cor_tissue$p.value,3)))
  tissue_PXS_GTEx_results <- tissue_PXS_GTEx_results %>% add_row(
    tissue = tissue, r = cor_tissue$estimate, p = cor_tissue$p.value
  )
}
tissue_PXS_order <- (tissue_PXS_GTEx_results %>% arrange(p))$tissue
ggplot(tissue_PXS_GTEx_results, aes(x = factor(tissue, levels = tissue_PXS_order), y = -log10(p))) +
  geom_col(aes(fill = p < 0.05, color = substring(tissue,1,5)=="Brain"), size=1.25) +
  geom_abline(slope = 0 , intercept = -log10(0.05)) +
  scale_color_manual(values=c("white", "black"),
                     breaks = c(FALSE, TRUE)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# looks at GTEx specific tissue data
for (i in 1:(nrow(jobID_tbl))) {
  jobID <- jobID_tbl$JobID[i]
  term <- jobID_tbl$term[i]
  traitname <- jobID_tbl$traitname[i]
  print(paste(i,jobID,term,traitname))
  
  loc_file <- paste0(dir_FUMA,"FUMA_job",jobID,"/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out")
  job_data <- as_tibble(fread(loc_file, skip="VARIABLE"))
  
  if (i==1) {GTEx_specific <- job_data %>% mutate(term = term, traitname = traitname)
  } else {
    GTEx_specific <- GTEx_specific %>%
      add_row(job_data %>% mutate(term = term, traitname = traitname))
  }
}
# Bonferroni correction (FUMA already Bonferonni corrects within each trait)
GTEx_specific$P_adjusted <- p.adjust(GTEx_specific$P / length(unique(GTEx_specific$VARIABLE)), method = 'bonferroni')
GTEx_specific$significant_P_adj <- GTEx_specific$P_adjusted < 0.05
GTEx_specific %>% filter(significant_P_adj) %>% arrange(P_adjusted) %>% print(n=Inf)
# tile plot
ggplot(GTEx_specific, aes(x = VARIABLE,
                         y = fct_rev(factor(traitname, levels=expo_order)))) +
  geom_tile(aes(fill=-log10(P_adjusted))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = length(expo_order)-0.5, ymax = length(expo_order)+0.5),
            color = "black", fill = NA) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=30), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=50), expand = c(0, 0)) +
  scale_fill_gradient2(low="white",mid="white", high="green", midpoint = -log10(0.05)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Tissue") + ylab("Trait") +
  labs(title = "Average gene expression per tissue for each behavior",
       subtitle = "Using MAGMA with GTEx through FUMA. White: p<0.05 (adjusted)")

#### GWAS CATALOG ####

# looks at GWAS_catalog data
for (i in 1:(nrow(jobID_tbl)-1)) {
  jobID <- jobID_tbl$JobID[i]
  term <- jobID_tbl$term[i]
  traitname <- jobID_tbl$traitname[i]
  print(paste(i,jobID,term,traitname))
  
  loc_file <- paste0(dir_FUMA,"FUMA_job",jobID,"/gwascatalog.txt")
  job_data <- as_tibble(fread(loc_file, skip="GenomicLocus"))
  if (nrow(job_data)==0) {next}
  job_data <- job_data %>% mutate(P = as.numeric(P))
  
  if (i==1) {GWAS_catalog <- job_data %>% mutate(term = term, traitname = traitname)
  } else {
    GWAS_catalog <- GWAS_catalog %>%
      add_row(job_data %>% mutate(term = term, traitname = traitname))
  }
}
GWAS_catalog %>% group_by(Trait) %>% summarize(n=n()) %>% arrange(-n)
