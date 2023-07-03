# Libraries and paths ####
library(GenomicRanges)
library(tidyverse)
library(data.table)
library(ggrepel)
library(rjson)
library(curl)
source('code/paths.R')

dir_FUMA <- paste0(dir_scratch,"FUMA_results/")
dir_figs <- paste0(dir_results,"figures/")

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

# gets order of exposures from decreasing p-value
coeffs <- as_tibble(fread(paste0(dir_script,"../input_data/PXS_coefficients.txt")))
expo_order <- factor(c(shortnames$shortname[1],
                       (coeffs %>%
                          filter(term %in% fields[fields$use_type=="exposure",]$term) %>% 
                          left_join(fields%>%select(term,shortname), by="term") %>%
                          arrange(p.value))$shortname) )

# shared theme for multiple plots
tileplot_theme <- theme(
  legend.key.size = unit(2, "mm"),
  legend.title = element_text(size=6),
  legend.text = element_text(size=5),
  legend.margin = margin(0,-1,0,-2, unit="mm"),
  legend.text.align = 1,
  plot.title = element_text(size=7),
  plot.subtitle = element_text(size=5),
  axis.title = element_text(size=6),
  axis.text = element_text(size=5)
)

# GTEx Data ####

## General Tissue data ####
for (i in 1:(nrow(jobID_tbl))) {
  jobID <- jobID_tbl$JobID[i]
  term <- jobID_tbl$term[i]
  shortname <- jobID_tbl$shortname[i]
  print(paste(i,jobID,term,shortname))
  
  loc_file <- paste0(dir_FUMA,"FUMA_job",jobID,"/magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out")
  job_data <- as_tibble(fread(loc_file, skip="VARIABLE"))
  
  if (i==1) {GTEx_general <- job_data %>% mutate(term = term, shortname = shortname)
  } else {
    GTEx_general <- GTEx_general %>%
      add_row(job_data %>% mutate(term = term, shortname = shortname))
  }
}
# Bonferroni correction (FUMA already Bonferonni corrects within each trait)
GTEx_general$P_adjusted <- p.adjust(GTEx_general$P / length(unique(GTEx_general$VARIABLE)), method = 'bonferroni')
GTEx_general$significant_P_adj <- GTEx_general$P_adjusted < 0.05
GTEx_general$tissue_name <- sapply(GTEx_general$VARIABLE, function(x) gsub("_", " ", x))
GTEx_general %>% filter(significant_P_adj) %>% arrange(P_adjusted) %>% print(n=Inf)
# saves to system
loc_out <- paste0(dir_results,"GTEx_general.txt")
fwrite(GTEx_general, loc_out, sep="\t")
### tile plot ####
ggplot(GTEx_general, aes(x = tissue_name,
                         y = fct_rev(factor(shortname, levels=expo_order)))) +
  geom_tile(aes(fill=-log10(P_adjusted))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = length(expo_order)-0.5, ymax = length(expo_order)+0.5),
            color = "black", fill = NA) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=30), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=20), expand = c(0, 0)) +
  scale_fill_gradient2(low="#DE1B1B",mid="white", high="#93CF1A", midpoint = -log10(0.05)) +
  xlab("Tissue Category") + ylab("Trait") +
  labs(#title = "Average gene expression per tissue category for each behavior",
       #subtitle = "Using MAGMA with GTEx through FUMA. White: adjusted p = 0.05",
       fill = bquote(-log[10](P)) ) +
  theme_light() +
  tileplot_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size=5))
# saves to system
loc_fig <- paste0(dir_figs,"GTEx_category")
ggsave(paste0(loc_fig,".png"), width=180, height=120, units="mm", dpi=300)
ggsave(paste0(loc_fig,".pdf"), width=180, height=120, units="mm", dpi=300)


## Specific tissue data ####
for (i in 1:(nrow(jobID_tbl))) {
  jobID <- jobID_tbl$JobID[i]
  term <- jobID_tbl$term[i]
  shortname <- jobID_tbl$shortname[i]
  print(paste(i,jobID,term,shortname))
  
  loc_file <- paste0(dir_FUMA,"FUMA_job",jobID,"/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out")
  job_data <- as_tibble(fread(loc_file, skip="VARIABLE"))
  
  if (i==1) {GTEx_specific <- job_data %>% mutate(term = term, shortname = shortname)
  } else {
    GTEx_specific <- GTEx_specific %>%
      add_row(job_data %>% mutate(term = term, shortname = shortname))
  }
}
# Bonferroni correction (FUMA already Bonferonni corrects within each trait)
GTEx_specific$P_adjusted <- p.adjust(GTEx_specific$P / length(unique(GTEx_specific$VARIABLE)), method = 'bonferroni')
GTEx_specific$significant_P_adj <- GTEx_specific$P_adjusted < 0.05
GTEx_specific$tissue_name <- sapply(GTEx_specific$FULL_NAME, function(x) gsub("_", " ", x))
GTEx_specific %>% filter(significant_P_adj) %>% arrange(P_adjusted) %>% print(n=Inf)
# saves to system
loc_out <- paste0(dir_results,"GTEx_specific.txt")
fwrite(GTEx_specific, loc_out, sep="\t")
### tile plot ####
ggplot(GTEx_specific, aes(x = tissue_name,
                         y = fct_rev(factor(shortname, levels=expo_order)))) +
  geom_tile(aes(fill=-log10(P_adjusted))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = length(expo_order)-0.5, ymax = length(expo_order)+0.5),
            color = "black", fill = NA) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=25), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=20), expand = c(0, 0)) +
  scale_fill_gradient2(low="#DE1B1B",mid="white", high="#93CF1A", midpoint = -log10(0.05)) +
  xlab("Tissue") + ylab("Trait") +
  labs(#title = "Average gene expression per tissue for each behavior",
       #subtitle = "Using MAGMA with GTEx through FUMA. White: adjusted p = 0.05",
       fill = bquote(-log[10](P)) ) +
  tileplot_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size=5))

# saves to system
loc_fig <- paste0(dir_figs,"GTEx_specific")
ggsave(paste0(loc_fig,".png"), width=180, height=120, units="mm", dpi=300)
ggsave(paste0(loc_fig,".pdf"), width=180, height=120, units="mm", dpi=300)

## Combined tissue data ####
# keeps list of tissues with at least one significant GTEx result
keep_general <- (GTEx_general %>% group_by(tissue_name) %>%
                   summarize(significant_P_adj = any(significant_P_adj)) %>%
                   filter(significant_P_adj))$tissue_name %>% unname()
keep_specific <- (GTEx_specific %>% group_by(tissue_name) %>%
                    summarize(significant_P_adj = any(significant_P_adj)) %>%
                    filter(significant_P_adj))$tissue_name %>% unname()
# keeps only the general tissue type if its duplicated as a specific tissue
keep_specific <-  keep_specific[!keep_specific %in% keep_general]

GTEx_combined <- GTEx_general %>% filter(tissue_name %in% keep_general) %>%
  add_row( GTEx_specific %>% select(-FULL_NAME) %>% filter(tissue_name %in% keep_specific) ) %>%
  select(-VARIABLE, -TYPE)

### tile plot ####
ggplot(GTEx_combined, aes(x = factor(tissue_name, levels = c(keep_general, keep_specific)),
                          y = fct_rev(factor(shortname, levels=expo_order)))) +
  geom_tile(aes(fill=-log10(P_adjusted))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = length(expo_order)-0.5, ymax = length(expo_order)+0.5),
            color = "black", fill = NA) +
  geom_rect(aes(xmin = 0.5, xmax = length(keep_general)+0.5),
                ymin = -Inf, ymax = Inf,
            color = "black", fill = NA) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=15), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=20), expand = c(0, 0)) +
  scale_fill_gradient2(low="#DE1B1B",mid="white", high="#93CF1A", midpoint = -log10(0.05)) +
  xlab("Tissue") + ylab("Trait") +
  labs(#title = "Average gene expression per tissue for each behavior",
       #subtitle = "Using MAGMA with GTEx through FUMA. White: adjusted p = 0.05",
       fill = bquote(-log[10](P))) +
  tileplot_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size=5))
# saves to system
loc_fig <- paste0(dir_figs,"GTEx_combined")
ggsave(paste0(loc_fig,".png"), width=140, height=120, units="mm", dpi=300)
ggsave(paste0(loc_fig,".pdf"), width=140, height=120, units="mm", dpi=300)












# GWAS CATALOG ####

# looks at GWAS_catalog data
for (i in 1:(nrow(jobID_tbl)-1)) {
  jobID <- jobID_tbl$JobID[i]
  term <- jobID_tbl$term[i]
  shortname <- jobID_tbl$shortname[i]
  print(paste(i,jobID,term,shortname))
  
  loc_file <- paste0(dir_FUMA,"FUMA_job",jobID,"/gwascatalog.txt")
  job_data <- as_tibble(fread(loc_file, skip="GenomicLocus"))
  if (nrow(job_data)==0) {next}
  job_data <- job_data %>% mutate(P = as.numeric(P))
  
  if (i==1) {GWAS_catalog <- job_data %>% mutate(term = term, shortname = shortname)
  } else {
    GWAS_catalog <- GWAS_catalog %>%
      add_row(job_data %>% mutate(term = term, shortname = shortname))
  }
}
GWAS_catalog_trait_count <- GWAS_catalog %>% group_by(Trait) %>% summarize(n=n()) %>% arrange(-n)
# # list of studies with many individual traits which represent similar concepts, removed
# GWAS_catalog %>% group_by(Study, Link) %>% summarize(n=n()) %>% arrange(-n)
# offender_studies <- c("www.ncbi.nlm.nih.gov/pubmed/35668104","www.ncbi.nlm.nih.gov/pubmed/34503513")

# #GWAS_catalog_traits <- GWAS_catalog_trait_count$Trait
# GWAS_catalog_trait_count_subset <- GWAS_catalog %>%
#   filter(!Link %in% offender_studies) %>%
#   group_by(Trait) %>% summarize(n = n()) %>% arrange(-n) %>%
#   filter(row_number() <= 150) #filter(n >= 2^6)
# GWAS_catalog_trait_count_subset$Trait %>% unique() %>% paste(collapse="|")
# # Fed this ^ into ChatGPT Plus (GPT4), and then manually fixed errors, check github repository

# goes through file of groupings
groupings_file <- readLines("input_data/GWAS_catalog_groupings.txt")
#groupings_file <- groupings_file[groupings_file != ""]
groupings_tbl <- tibble(Group = as.character(), Trait = as.character())
for (line in groupings_file) {
  Group <- strsplit(line, ": ")[[1]][1]
  traits_text <- strsplit(line, ": ")[[1]][2]
  traits <- strsplit(traits_text,", ")[[1]]
  groupings_tbl <- groupings_tbl %>% add_row(Group = Group, Trait = traits)
}
# manually fixes cases where a ", " is in the trait name
groupings_tbl$Trait <- gsub("\\\\", "", groupings_tbl$Trait)
# groupings_tbl$Group %>% unique()
# # check for duplicates
# groupings_tbl %>% 
#   filter(Trait %in% (groupings_tbl%>%group_by(Trait)%>%summarize(n=n()) %>% filter(n>1))$Trait) %>%
#   arrange(Trait)
# # check for missingness
# (missing_traits <- (GWAS_catalog_trait_count_subset %>% filter(!Trait %in% groupings_tbl$Trait))$Trait)

# saves groupings_tbl to system
loc_out <- paste0(dir_results,"GWAS_catalog_groupings.txt")
fwrite(groupings_tbl, loc_out, sep="\t")

GWAS_catalog <- GWAS_catalog %>% left_join(groupings_tbl, by="Trait")
GWAS_catalog_group_count <- GWAS_catalog %>%
  group_by(Group) %>% summarize(n = n()) %>%
  filter(!is.na(Group)) %>% arrange(-n)

# saves tables
loc_out <- paste0(dir_results,"GWAS_catalog_group_count.txt")
fwrite(GWAS_catalog_group_count, loc_out, sep="\t")

loc_out <- paste0(dir_results,"GWAS_catalog_trait_count.txt")
fwrite(GWAS_catalog_trait_count, loc_out, sep="\t")




# Genes ####

## extracts FUMA gene results ####
for (i in 1:(nrow(jobID_tbl))) {
  jobID <- jobID_tbl$JobID[i]
  term <- jobID_tbl$term[i]
  shortname <- jobID_tbl$shortname[i]
  print(paste(i,jobID,term,shortname))
  
  loc_file <- paste0(dir_FUMA,"FUMA_job",jobID,"/genes.txt")
  job_data <- as_tibble(fread(loc_file, skip="ensg"))
  
  if (i==1) {Genes <- job_data %>% mutate(term = term, shortname = shortname,
                                          GenomicLocus = as.character(GenomicLocus))
  } else {
    Genes <- Genes %>%
      add_row(job_data %>% mutate(term = term, shortname = shortname,
                                  GenomicLocus = as.character(GenomicLocus)))
  }
}
# gets T2D gene associations from AMP
AMP_T2D_genes <- as_tibble(fread("input_data/AMP_T2D_gene_table.csv"))
# makes joint list of gene associations
Bonferroni_gene <- 2.5E-6 # 0.05 / 20000
Genes_PXS <- Genes %>% filter(term=="PXS_T2D") %>%
  select(gene = symbol, p_PXS_T2D = minGwasP) %>%
  unique() %>%
  full_join(AMP_T2D_genes %>% select(gene, p_T2D = minP), by="gene") %>%
  mutate(p_T2D = ifelse(is.na(p_T2D),1,p_T2D),
         p_PXS_T2D = ifelse(is.na(p_PXS_T2D),1,p_PXS_T2D),
         sig_T2D = p_T2D < Bonferroni_gene,
         sig_PXS_T2D = p_PXS_T2D < Bonferroni_gene)
# compares significance
(gene_count <- table(Genes_PXS %>% select(starts_with("sig_"))) )
gene_count[2,2] / sum(gene_count[,2])

# gets table of most repreated genes
Genes %>% group_by(symbol) %>% summarize(n=n()) %>% arrange(-n)
n6_Genes <- (Genes %>% group_by(symbol) %>% summarize(n=n()) %>%
             filter(n>5) %>% arrange(-n))$symbol
Genes_top <- Genes %>% select(symbol,term,shortname) %>%
  add_row(AMP_T2D_genes %>% filter(minP < Bonferroni_gene) %>% select(symbol = gene) %>%
            mutate(term = "T2D", shortname = "Type 2 Diabetes")) %>%
  filter(symbol %in% n6_Genes) %>%
  mutate(gene = as.factor(symbol))

# plots crappy heatmap of repeated genes
ggplot(Genes_top, aes(x = factor(gene, levels = n6_Genes),
                       y = fct_rev(factor(shortname, levels=c("Type 2 Diabetes", levels(expo_order)[expo_order])) ))) +
  geom_tile(fill = "red") + theme_light()

# Number of genes by trait
Genes_trait <- Genes %>% filter(minGwasP < Bonferroni_gene) %>%
  group_by(shortname) %>% summarize(n = n()) %>% arrange(-n)
# saves to system
loc_out <- paste0(dir_results,"gene_associations_n.txt")
fwrite(Genes_trait, loc_out, sep="\t")

# Genomic Loci ####
for (i in 1:(nrow(jobID_tbl))) {
  jobID <- jobID_tbl$JobID[i]
  term <- jobID_tbl$term[i]
  shortname <- jobID_tbl$shortname[i]
  print(paste(i,jobID,term,shortname))
  
  loc_file <- paste0(dir_FUMA,"FUMA_job",jobID,"/GenomicRiskLoci.txt")
  job_data <- as_tibble(fread(loc_file, skip="GenomicLocus"))
  
  if (i==1) {GRLoci <- job_data %>% mutate(term = term, shortname = shortname)
  } else {
    GRLoci <- GRLoci %>%
      add_row(job_data %>% mutate(term = term, shortname = shortname))
  }
}
# shared loci
gr <- GenomicRanges::GRanges(seqnames = GRLoci$chr,
                             ranges = IRanges(start = GRLoci$start, end = GRLoci$end)) %>%
  GenomicRanges::reduce()
GRLoci_merged <- as_tibble(gr) %>%
  mutate(GLM_ID = row_number()) %>%
  select(GLM_chr = seqnames, GLM_start = start, GLM_end = end, GLM_ID)

GRLoci <- GRLoci %>% arrange(chr, start) %>%
  mutate(GLM_ID=0)
for (i in 1:nrow(GRLoci)) {
  GRLoci_merged_chr <- GRLoci_merged %>% filter(GLM_chr==GRLoci$chr[i])
  for (ii in 1:nrow(GRLoci_merged_chr)) {
    if (GRLoci$start[i] >= GRLoci_merged_chr$GLM_start[ii] & GRLoci$end[i] <= GRLoci_merged_chr$GLM_end[ii]) {
      GRLoci$GLM_ID[i] <- GRLoci_merged_chr$GLM_ID[ii]
      break
    }
  }
  if (i %% 50 == 0) {print(i)}
}
GRLoci <- GRLoci %>% left_join(GRLoci_merged, by="GLM_ID")
GRLoci %>% group_by(GLM_ID, GLM_chr, GLM_start, GLM_end) %>% summarize(n=n()) %>% arrange(-n)
n3_GLM <- (GRLoci %>% group_by(GLM_ID) %>% summarize(n=n()) %>%
             filter(n>2) %>% arrange(-n))$GLM_ID
GRLoci_top <- GRLoci %>% filter(GLM_ID %in% n3_GLM) %>%
  mutate(GLM_ID = as.factor(GLM_ID))

# plots crappy heatmap of repeated loci
ggplot(GRLoci_top, aes(x = factor(GLM_ID, levels = n3_GLM),
                       y = fct_rev(factor(shortname, levels=expo_order)))) +
  geom_tile(fill = "red") + theme_light()

# Number of loci by traits
GRLoci %>% group_by(shortname) %>% summarize(n = n()) %>% arrange(-n)
# prints table
loc_out <- paste0(dir_results,"genomic_risk_loci.txt")
fwrite(GRLoci, loc_out, sep="\t")
