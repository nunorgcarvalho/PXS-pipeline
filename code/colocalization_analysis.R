library(tidyverse)
library(data.table)
library(coloc)

# Shared Code ####
source('code/paths.R')

dir_FUMA <- paste0(dir_scratch,"FUMA_results/")
dir_coloc <- paste0(dir_results,"figures/colocs/")
dir.create(dir_coloc, showWarnings = FALSE)
BP_window <- 50000 # used both above and used below

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

# function for transforming large LMM summary file to coloc dataset within one region
convert2dataset <- function(LMM, chromosome, BP_min, BP_max, N=300000, type="quant") {
  # filters to specified region and defines variance (~ SE^2)
  LMM_region <- LMM %>% filter(CHR==chromosome, between(BP, BP_min, BP_max)) %>%
    mutate(varbeta = SE^2)
  
  # removes SNPs with duplicate results (due to multiple alleles)
  dupe_SNPs <- (LMM_region %>% group_by(SNP) %>% summarize(n=n()) %>% filter(n>1))$SNP
  LMM_region <- LMM_region %>% filter(!SNP %in% dupe_SNPs)
  
  # formats list correctly for coloc
  LMM_dataset <- list(
    snp = LMM_region$SNP,
    position = LMM_region$BP,
    MAF = LMM_region$A1FREQ,
    beta = LMM_region$BETA,
    varbeta = LMM_region$varbeta,
    N = N,
    type = type
  )
  return(LMM_dataset)
}

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
Genes %>% select(symbol, term) %>% distinct() %>% 
  group_by(symbol) %>% summarize(n=n()) %>% arrange(-n)
# extracts genes found significant in at least 7 traits
genes_of_interest <- ( Genes %>% select(symbol, term) %>% distinct() %>% 
  group_by(symbol) %>% summarize(n=n()) %>% filter(n>=7))$symbol


coloc_datasets <- vector("list",length = nrow(shortnames))
# loop through all terms ####
for (i in 1:nrow(shortnames)) {
  term <- shortnames$term[i]
  # sets summary file path
  if (term == "PXS_T2D") {loc_LMM <- "scratch/T2D/LMM_PXS_T2D_bgen.txt"
  } else {loc_LMM <- paste0("scratch/exposures/",term,"/LMM_",term,"_bgen.txt")}
  # reads summary file
  print(paste0(i," :: Reading: ", term))
  LMM <- as_tibble(fread(loc_LMM))
  
  coloc_datasets[[i]] <- vector("list",length = length(genes_of_interest))
  # loops through each gene
  for (ii in 1:length(genes_of_interest)) {
    # gets gene location/window
    gene <- genes_of_interest[ii]
    chromosome <- median((Genes %>% filter(symbol == gene))$chr)
    BP_min <- min((Genes %>% filter(symbol == gene))$start) - BP_window
    BP_max <- max((Genes %>% filter(symbol == gene))$end) + BP_window
    # converts/filters LMM to gene region
    print(paste0(i, " :: ", gene, " :: Converting to dataset: ", term))
    LMM_dataset <- convert2dataset(LMM, chromosome, BP_min, BP_max)
    # adds filtered LMM to running nested list
    coloc_datasets[[i]][[ii]] <- LMM_dataset
  }
}
# saves filtered nested list to system (saves time if running again)
#saveRDS(coloc_datasets, file = "scratch/coloc_datasets.rds")

# colocalization analysis ####
coloc_tibble <- tibble(gene = as.character(),
                       term1 = as.character(),
                       term2 = as.character(),
                       H4PP = as.numeric())
gene_tibble <- tibble(gene = as.character(),
                      term = as.character(),
                      sig_gene = as.logical())
# loops throug each gene
for (k in 1:length(genes_of_interest)) {
  gene <- genes_of_interest[k]
  # loops through each unique combination of 2 traits
  for (i in 1:nrow(shortnames)) {
    for (j in 1:i) {
      if (i == j) {next}
      print(paste(i,j))
      # runs colocalization analysis using bayes factor
      coloc1 <- coloc.abf(coloc_datasets[[i]][[k]], coloc_datasets[[j]][[k]])
      # saves H4 posterior probability
      coloc_tibble <- coloc_tibble %>%
        add_row(
          gene = gene,
          term1 = shortnames$term[i],
          term2 = shortnames$term[j],
          H4PP = coloc1$summary[6]
        )
    }
    # checks if gene is considered significant for this trait
    if (nrow(Genes %>% filter(symbol == gene, term == shortnames$term[i]) > 0 )) {
      sig_gene <- TRUE
    } else {sig_gene <- FALSE}
    # adds to running table
    gene_tibble <- gene_tibble %>%
      add_row(gene = gene, term = shortnames$term[i], sig_gene = sig_gene)
  }
}

# names terms after their shortened names
coloc_tibble <- coloc_tibble %>%
  mutate(term1 = factor(term1, levels = shortnames$term, labels = shortnames$shortname)) %>%
  mutate(term2 = factor(term2, levels = shortnames$term, labels = shortnames$shortname))

# tileplots ####
# loops through each gene
for (gene in genes_of_interest) {
  # gets vector denoting which traits to bold in axes
  boldings <- ( gene_tibble %>% filter(gene == .GlobalEnv$gene) )$sig_gene
  # makes  plot
  gg <- ggplot(coloc_tibble %>% filter(gene == .GlobalEnv$gene),
         aes(x=term1, y = term2)) +
    geom_tile(aes(fill = H4PP)) +
    geom_text(aes(label = sprintf("%.2f", H4PP)), color="black", size=2) +
    scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=20), expand = c(0, 0)) +
    scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=20), expand = c(0, 0)) +
    scale_fill_gradientn(colors=c("#DE1B1B","white", "#93CF1A"), values = c(0,0.5,1)) +
    theme_light() +
    theme(axis.text = element_text(size=5),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,
                                 face = ifelse(boldings[-1], "bold", "plain")),
      axis.text.y = element_text(face = ifelse(boldings, "bold", "plain")),
      legend.key.size = unit(2, "mm"),
      legend.title = element_text(size=6),
      legend.text = element_text(size=5),
      legend.margin = margin(0,-1,0,-2, unit="mm"),
      legend.text.align = 1,
      plot.title = element_text(size=7),
      plot.subtitle = element_text(size=5),
      axis.title = element_text(size=6),
    ) +
    labs(title = paste0("Posterior probability of shared causal variant in ",gene," gene region") ,
         subtitle = "Axis text that is bolded denotes traits where gene is considered significantly associated (p < 2.5E-6)",
         fill = "Prob.") +
    xlab("Trait") + ylab("Trait")
  
  # saves plot to system
  loc_fig <- paste0(dir_coloc, "coloc_",gene)
  ggsave(paste0(loc_fig,".png"), width=180, height=120, units="mm", dpi=300)
  #ggsave(paste0(loc_fig,".pdf"), width=180, height=120, units="mm", dpi=300)
}
