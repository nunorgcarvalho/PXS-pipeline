# Libraries and paths ####
library(GenomicRanges)
library(tidyverse)
library(data.table)
library(ggrepel)
library(rjson)
library(curl)
source('code/00_paths.R')

dir_FUMA <- paste0(dir_scratch,"FUMA_results/")
dir.create(paste0(dir_results,"figures/"), showWarnings = FALSE)
dir_figs <- paste0(dir_results,"figures/")

# Shared code ####

## unzips and renames folder ####
zips <- list.files(paste0(dir_FUMA,'zips/'))
dir_FUMAs <- c()
for (i in 1:length(zips)) {
  zip <- zips[i]
  print(paste(i, 'Unzipping', zip))
  zip_path <- paste0(dir_FUMA,'zips/',zip)
  # extracts variable ran thru FUMA
  con <- unz(zip_path, "params.config")
  params <- readLines(con)
  close(con)
  term <- params[str_detect(params, 'title = ')] %>%
    str_replace_all('title = ', '') %>% # assumes term is in title
    str_replace_all('LMM.ALL.', '') # for later consistency
  dir_FUMAs <- c(dir_FUMAs, term)
  
  # saves to new folder unless it already exists
  dir_term <- paste0(dir_FUMA,'FUMA_',term)
  if (!dir.exists(dir_term)) {
    unzip(zip_path, exdir = dir_term)
  }
}


# used for plotting purposes
sumtbl <- as_tibble(fread(paste0(dir_scratch,'general_results/sumtbl.tsv')))
bvr_order <- factor( (sumtbl %>% arrange(-abs(beta_norm)))$shortname )
term_BRS <- sumtbl$term[1]

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


# function for reading FUMA results for all jobs
read_FUMAs <- function(dir_FUMAs, FUMA_file, skip='VARIABLE') {
  for (i in 1:length(dir_FUMAs)) {
    term <- dir_FUMAs[i]
    print(paste(i,term))
    
    loc_file <- paste0(dir_FUMA,'FUMA_',term,'/',FUMA_file)
    if (!file.exists(loc_file)) { # skips missing file
      print(paste(i, term, FUMA_file, 'file not fond'))
      next
    }
    job_data <- as_tibble(fread(loc_file, skip=skip))
    if (nrow(job_data)==0) { # skips empty data table
      print(paste(i, term, FUMA_file, 'data table is empty'))
      next
    }
    if (FUMA_file == 'genes.txt') { # needed for X,Y chromosomes
      job_data <- job_data %>% mutate(GenomicLocus = as.character(GenomicLocus))
    }
    
    if (i==1) {FUMA_data <- job_data %>% mutate(term = term)
    } else {
      FUMA_data <- FUMA_data %>% add_row(job_data %>% mutate(term = term))
    }
  }
  return(FUMA_data)
}

# GTEx Data ####

## GTEx common theme ####
GTEx_theme <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                    legend.key.size = unit(4, "mm"),
                    legend.title = element_text(size=5))

## General Tissue data ####
GTEx_general <- read_FUMAs(dir_FUMAs,
                           'magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out') %>%
  # restricts to just individual behaviors and BRS
  filter(term %in% sumtbl$term) %>%
  left_join(sumtbl %>% select(term,shortname),by='term') %>%
  mutate(tissue_name = sapply(VARIABLE, function(x) gsub("_", " ", x)))

# Bonferroni correction (FUMA already Bonferonni corrects within each trait)
GTEx_general$P_adjusted <- p.adjust(GTEx_general$P / length(unique(GTEx_general$VARIABLE)), method = 'bonferroni')
GTEx_general$significant_P_adj <- GTEx_general$P_adjusted < 0.05
GTEx_general %>% filter(significant_P_adj) %>% arrange(P_adjusted) %>%
  select(tissue_name, term, shortname,P_adjusted) %>% print(n=Inf)
# saves to system
loc_out <- paste0(dir_scratch, 'general_results/GTEx_general.tsv')
fwrite(GTEx_general, loc_out, sep="\t")
### tile plot ####
ggplot(GTEx_general, aes(x = tissue_name,
                         y = fct_rev(factor(shortname, levels=bvr_order)))) +
  geom_tile(aes(fill=-log10(P_adjusted))) +
  # geom_rect(aes(xmin = -Inf, xmax = Inf,
  #               ymin = length(bvr_order)-0.5, ymax = length(bvr_order)+0.5),
  #           color = "black", fill = NA) +
  scale_x_discrete(labels = function(x) str_wrap(x, width=30), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width=20), expand = c(0, 0)) +
  scale_fill_gradient2(low="#DE1B1B",mid="white", high="#93CF1A", midpoint = -log10(0.05)) +
  xlab("Tissue Category") + ylab("Trait") +
  labs(#title = "Average gene expression per tissue category for each behavior",
       #subtitle = "Using MAGMA with GTEx through FUMA. White: adjusted p = 0.05",
       fill = bquote(-log[10](P)) ) +
  theme_light() +
  tileplot_theme +
  GTEx_theme
# saves to system
loc_fig <- paste0(dir_figs,"GTEx_category")
ggsave(paste0(loc_fig,".png"), width=180, height=120, units="mm", dpi=300)
ggsave(paste0(loc_fig,".pdf"), width=180, height=120, units="mm", dpi=300)




## Specific tissue data ####
GTEx_specific <- read_FUMAs(dir_FUMAs,
                             'magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out') %>%
  # restricts to just individual behaviors and BRS
  filter(term %in% sumtbl$term) %>%
  left_join(sumtbl %>% select(term,shortname),by='term') %>%
  mutate(tissue_name = sapply(VARIABLE, function(x) gsub("_", " ", x)))

# Bonferroni correction (FUMA already Bonferonni corrects within each trait)
GTEx_specific$P_adjusted <- p.adjust(GTEx_specific$P / length(unique(GTEx_specific$VARIABLE)), method = 'bonferroni')
GTEx_specific$significant_P_adj <- GTEx_specific$P_adjusted < 0.05
GTEx_specific$tissue_name <- sapply(GTEx_specific$FULL_NAME, function(x) gsub("_", " ", x))
GTEx_specific %>% filter(significant_P_adj) %>% arrange(P_adjusted) %>%
  select(tissue_name, term, shortname,P_adjusted) %>% print(n=Inf)
# saves to system
loc_out <- paste0(dir_scratch, 'general_results/GTEx_specific.tsv')
fwrite(GTEx_specific, loc_out, sep="\t")
### tile plot ####
ggplot(GTEx_specific, aes(x = tissue_name,
                         y = fct_rev(factor(shortname, levels=bvr_order)))) +
  geom_tile(aes(fill=-log10(P_adjusted))) +
  # geom_rect(aes(xmin = -Inf, xmax = Inf,
  #               ymin = length(bvr_order)-0.5, ymax = length(bvr_order)+0.5),
  #           color = "black", fill = NA) +
  scale_x_discrete(labels = function(x) str_wrap(x, width=25), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width=20), expand = c(0, 0)) +
  scale_fill_gradient2(low="#DE1B1B",mid="white", high="#93CF1A", midpoint = -log10(0.05)) +
  xlab("Tissue") + ylab("Trait") +
  labs(#title = "Average gene expression per tissue for each behavior",
       #subtitle = "Using MAGMA with GTEx through FUMA. White: adjusted p = 0.05",
       fill = bquote(-log[10](P)) ) +
  theme_light() +
  tileplot_theme +
  GTEx_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4))

# saves to system
loc_fig <- paste0(dir_figs,"GTEx_specific")
ggsave(paste0(loc_fig,".png"), width=180, height=120, units="mm", dpi=300)
ggsave(paste0(loc_fig,".pdf"), width=180, height=120, units="mm", dpi=300)

## Combined tissue data ####
# keeps list of tissues with at least one significant GTEx result for BRS
keep_general <- (GTEx_general %>% filter(term==term_BRS) %>% filter(significant_P_adj))$tissue_name %>% unname()
keep_specific <- (GTEx_specific %>% filter(term==term_BRS) %>% filter(significant_P_adj))$tissue_name %>% unname()
# keeps only the general tissue type if its duplicated as a specific tissue
keep_specific <-  keep_specific[!keep_specific %in% keep_general]

GTEx_combined <- GTEx_general %>% filter(tissue_name %in% keep_general) %>%
  add_row( GTEx_specific %>% select(-FULL_NAME) %>% filter(tissue_name %in% keep_specific) ) %>%
  select(-VARIABLE, -TYPE)

### tile plot ####
ggplot(GTEx_combined, aes(x = factor(tissue_name, levels = c(keep_general, keep_specific)),
                          y = fct_rev(factor(shortname, levels=bvr_order)))) +
  geom_tile(aes(fill=-log10(P_adjusted))) +
  # geom_rect(aes(xmin = -Inf, xmax = Inf,
  #               ymin = length(bvr_order)-0.5, ymax = length(bvr_order)+0.5),
  #           color = "black", fill = NA) +
  # geom_rect(aes(xmin = 0.5, xmax = length(keep_general)+0.5),
  #               ymin = -Inf, ymax = Inf,
  #           color = "black", fill = NA) +
  geom_segment(x = length(keep_general)+0.5, xend=length(keep_general)+0.5,
               y = -Inf, yend = Inf, color = "white", size=2) +
  geom_rect(aes(xmin = 0.5, xmax = length(keep_general)+0.4),
            ymin = -Inf, ymax = Inf,
            color = "black", fill = NA, size=0.3) +
  geom_rect(aes(xmin = length(keep_general)+0.6, xmax=Inf),
            ymin = -Inf, ymax = Inf,
            color = "black", fill = NA, size=0.3) +
  scale_x_discrete(labels = function(x) str_wrap(x, width=15), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width=20), expand = c(0, 0)) +
  scale_fill_gradient2(low="#DE1B1B",mid="white", high="#93CF1A", midpoint = -log10(0.05)) +
  xlab("Tissue") + ylab("Trait") +
  labs(#title = "Average gene expression per tissue for each behavior",
       #subtitle = "Using MAGMA with GTEx through FUMA. White: adjusted p = 0.05",
       fill = bquote(-log[10](P))) +
  tileplot_theme +
  GTEx_theme
# saves to system
loc_fig <- paste0(dir_figs,"GTEx_combined")
ggsave(paste0(loc_fig,".png"), width=140, height=120, units="mm", dpi=300)
ggsave(paste0(loc_fig,".pdf"), width=140, height=120, units="mm", dpi=300)

# Genes ####

## extracts FUMA gene results ####
Genes <- read_FUMAs(dir_FUMAs, 'genes.txt',skip='ensg') %>%
  # restricts to just individual behaviors and BRS
  filter(term %in% sumtbl$term) %>%
  left_join(sumtbl %>% select(term,shortname),by='term')

# gets T2D gene associations from AMP
AMP_T2D_genes <- as_tibble(fread("scratch_2023/AMP/AMP_T2D_gene_table.csv"))
# makes joint list of gene associations
Bonferroni_gene <- 2.5E-6 # 0.05 / 20000
Genes_BRS <- Genes %>% filter(term==term_BRS) %>%
  select(gene = symbol, p_BRS = minGwasP) %>%
  unique() %>%
  full_join(AMP_T2D_genes %>% select(gene, p_T2D = minP), by="gene") %>%
  mutate(p_T2D = ifelse(is.na(p_T2D),1,p_T2D),
         p_BRS = ifelse(is.na(p_BRS),1,p_BRS),
         sig_T2D = p_T2D < Bonferroni_gene,
         sig_BRS = p_BRS < Bonferroni_gene)
# compares significance
(gene_count <- table(Genes_BRS %>% select(starts_with("sig_"))) )
gene_count[2,2] / sum(gene_count[,2])

# gets table of most repeated genes
(Genes %>% group_by(symbol) %>% summarize(n=n()) %>% arrange(-n))$n %>% table()
n10_Genes <- (Genes %>% group_by(symbol) %>% summarize(n=n()) %>%
             filter(n>9) %>% arrange(-n))$symbol
Genes_top <- Genes %>% select(symbol,term,shortname) %>%
  add_row(AMP_T2D_genes %>%
            filter(minP < Bonferroni_gene) %>%
            select(symbol = gene) %>%
            mutate(term = "T2D",
                   shortname = "Type 2 Diabetes")) %>%
  filter(symbol %in% n10_Genes) %>%
  mutate(gene = as.factor(symbol))

# plots crappy heatmap of repeated genes
ggplot(Genes_top, aes(x = factor(gene, levels = n10_Genes),
                       y = fct_rev(factor(shortname, levels=c("Type 2 Diabetes", levels(bvr_order)[bvr_order])) ))) +
  geom_tile(fill = "red") + theme_light()

# Number of genes by trait
Genes_trait <- Genes %>% filter(minGwasP < Bonferroni_gene) %>%
  group_by(shortname) %>% summarize(n = n()) %>% arrange(-n)
# saves to system
loc_out <- paste0(dir_scratch, 'general_results/gene_associations_n.tsv')
fwrite(Genes_trait, loc_out, sep="\t")



# Genomic Loci ####
GRLoci <- read_FUMAs(dir_FUMAs, 'GenomicRiskLoci.txt', skip='GenomicLocus') %>%
  # restricts to just individual behaviors and BRS
  filter(term %in% sumtbl$term) %>%
  left_join(sumtbl %>% select(term,shortname),by='term')

# shared loci
gr <- GenomicRanges::GRanges(seqnames = GRLoci$chr,
                             ranges = IRanges(start = GRLoci$start, end = GRLoci$end)) %>%
  GenomicRanges::reduce()
GRLoci_merged <- as_tibble(gr) %>%
  mutate(GLM_ID = row_number()) %>%
  select(GLM_chr = seqnames, GLM_start = start, GLM_end = end, GLM_ID)

GRLoci <- GRLoci %>% arrange(chr, start) %>% mutate(GLM_ID=0)
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
(GRLoci %>% group_by(GLM_ID, GLM_chr, GLM_start, GLM_end) %>% summarize(n=n()) %>% arrange(-n))$n %>% table()
n6_GLM <- (GRLoci %>% group_by(GLM_ID) %>% summarize(n=n()) %>%
             filter(n>5) %>% arrange(-n))$GLM_ID
GRLoci_top <- GRLoci %>% filter(GLM_ID %in% n6_GLM) %>%
  mutate(GLM_ID = as.factor(GLM_ID))

# plots crappy heatmap of repeated loci
ggplot(GRLoci_top, aes(x = factor(GLM_ID, levels = n6_GLM),
                       y = fct_rev(factor(shortname, levels=bvr_order)))) +
  geom_tile(fill = "red") + theme_light()

# Number of loci by traits
GRLoci %>% group_by(shortname) %>% summarize(n = n()) %>% arrange(-n)
# prints table
loc_out <- paste0(dir_scratch, 'general_results/genomic_risk_loci.tsv')
fwrite(GRLoci, loc_out, sep="\t")




## BRS Loci ####
# identifies BRS loci not found in other traits
GLM_BRS <- GRLoci %>% filter(term==term_BRS)
GLM_BRS2 <- GRLoci %>% group_by(GLM_ID) %>%
  summarize(n_behaviors_sig = n() - 1) %>%
  filter(GLM_ID %in% GLM_BRS$GLM_ID) %>%
  left_join(GLM_BRS %>% select(GenomicLocus, GLM_ID), by="GLM_ID")


GL_BRS <- Genes %>% filter(term == term_BRS) %>%
  mutate(GenomicLocus = as.numeric(GenomicLocus)) %>%
  group_by(GenomicLocus) %>%
  summarize(n_genes=n(), minGwasP = min(minGwasP),
            genes = paste(unique(symbol), collapse = ", ")) %>%
  left_join(GRLoci %>% filter(term == term_BRS) %>% select(GenomicLocus, rsID, chr, pos, p), by="GenomicLocus") %>%
  mutate(concordant = (minGwasP == p)) %>% arrange(p) %>%
  select(-chr, -pos)

# contains allele, MAF, and beta data
SNPs <- as_tibble(fread(paste0(dir_FUMA, 'FUMA_', term_BRS,'/snps.txt'))) %>%
  select(rsID, chr, pos, non_effect_allele, effect_allele, MAF, beta, se, nearestGene) %>%
  mutate(nearestGene = str_replace_all(nearestGene, ':',', '))

GL_BRS2 <- GL_BRS %>%
  left_join(SNPs, by="rsID") %>%
  left_join(GLM_BRS2, by="GenomicLocus") %>%
  select(SNP=rsID, CHR = chr, BP = pos, A0 = non_effect_allele, A1 = effect_allele,
         MAF, beta, beta_se = se, p, closest_genes = genes, n_behaviors_sig)


loc_out <- paste0(dir_scratch, 'general_results/GenomicLoci_BRS.tsv')
fwrite(GL_BRS2, loc_out, sep="\t")

## BRS loci2 ####

# contains top SNPs
top_SNPs <- as_tibble(fread(paste0(dir_FUMA, 'FUMA_', term_BRS,'/IndSigSNPs.txt'))) %>%
  group_by(GenomicLocus) %>%
  filter(p == min(p)) %>% arrange(p) %>%
  left_join(SNPs, by=c('rsID','chr','pos'))

BRS_Genes <- as_tibble(fread(paste0(dir_FUMA, 'FUMA_', term_BRS,'/genes.txt'))) %>%
  select(GenomicLocus, symbol, minGwasP) %>%
  add_row(
    top_SNPs %>%
      select(GenomicLocus, symbol = nearestGene, minGwasP = p) %>%
      ungroup()
  ) %>%
  group_by(GenomicLocus) %>%
  arrange(minGwasP, symbol) %>%
  unique() %>%
  filter(!grepl(', ', symbol)) %>%
  summarize(genes = paste(unique(symbol), collapse = ", "))

BRS_tableS3 <- top_SNPs %>%
  ungroup() %>%
  left_join(BRS_Genes, by='GenomicLocus') %>%
  select(topSNP = rsID, CHR=chr, BP=pos, A0=non_effect_allele, A1=effect_allele,
         MAF, beta, se, p, closest_gene = nearestGene, all_locus_genes = genes)

loc_out <- paste0(dir_scratch, 'general_results/GenomicLoci_BRS.csv')
fwrite(BRS_tableS3, loc_out)
