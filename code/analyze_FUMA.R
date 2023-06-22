# Libraries and paths ####
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

# GTEx Data ####
# gets order of exposures from decreasing p-value
coeffs <- as_tibble(fread(paste0(dir_script,"../input_data/PXS_coefficients.txt")))
expo_order <- factor(c("PXS for Type 2 Diabetes onset",
                       (coeffs %>%
                          filter(term %in% fields[fields$use_type=="exposure",]$term) %>% 
                          left_join(fields%>%select(term,traitname), by="term") %>%
                          arrange(p.value))$traitname) )

## GTEx Category data ####
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
### tile plot ####
ggplot(GTEx_general, aes(x = VARIABLE,
                         y = fct_rev(factor(traitname, levels=expo_order)))) +
  geom_tile(aes(fill=-log10(P_adjusted))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = length(expo_order)-0.5, ymax = length(expo_order)+0.5),
            color = "black", fill = NA) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=30), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=50), expand = c(0, 0)) +
  scale_fill_gradient2(low="red",mid="white", high="green", midpoint = -log10(0.05)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Tissue Category") + ylab("Trait") +
  labs(title = "Average gene expression per tissue category for each behavior",
       subtitle = "Using MAGMA with GTEx through FUMA. White: p<0.05 (Bonferroni adjusted)")

## GTEx Cpecific tissue data ####
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
### tile plot ####
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

# GWAS CATALOG ####

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
# list of studies with many individual traits which represent similar concepts, removed
GWAS_catalog %>% group_by(Study, Link) %>% summarize(n=n()) %>% arrange(-n)
offender_studies <- c("www.ncbi.nlm.nih.gov/pubmed/35668104","www.ncbi.nlm.nih.gov/pubmed/34503513")
GWAS_catalog_trait_count <- GWAS_catalog %>% group_by(Trait) %>% summarize(n=n()) %>% arrange(-n)
#GWAS_catalog_traits <- GWAS_catalog_trait_count$Trait
GWAS_catalog_trait_count_subset <- GWAS_catalog %>%
  filter(!Link %in% offender_studies) %>%
  group_by(Trait) %>% summarize(n = n()) %>% arrange(-n) %>%
  filter(n >= 2^6)
GWAS_catalog_trait_count_subset$Trait %>% unique() %>% paste(collapse="|")
# Fed this ^ into ChatGPT Plus (GPT4)

# goes through file of groupings
groupings_file <- readLines("../input_data/GWAS_catalog_groupings.txt")
groupings_file <- groupings_file[groupings_file != ""]
groupings_tbl <- tibble(Group = as.character(), Trait = as.character())
for (line in groupings_file) {
  Group <- strsplit(line, ": ")[[1]][1]
  traits_text <- strsplit(line, ": ")[[1]][2]
  traits <- strsplit(traits_text,", ")[[1]]
  groupings_tbl <- groupings_tbl %>% add_row(Group = Group, Trait = traits)
}
groupings_tbl$Group %>% unique()
# check for duplicates
groupings_tbl %>% 
  filter(Trait %in% (groupings_tbl%>%group_by(Trait)%>%summarize(n=n()) %>% filter(n>1))$Trait) %>%
  arrange(Trait)
# check for missingness
(missing_traits <- (GWAS_catalog_trait_count_subset %>% filter(!Trait %in% groupings_tbl$Trait))$Trait)
# manually fixes groupings
groupings_tbl <- groupings_tbl[-c(45,66,84,134),]
groupings_tbl %>%
  add_row(Group = c("Metabolic Traits","Cardiovascular Traits","Cardiovascular Traits"),
          Trait = missing_traits)

GWAS_catalog <- GWAS_catalog %>% left_join(groupings_tbl, by="Trait")
GWAS_catalog %>%
  group_by(Group) %>% summarize(n = n()) %>%
  filter(!is.na(Group)) %>% arrange(-n)
GWAS_catalog_trait_count
