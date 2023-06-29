topGenes <- Genes %>%
  group_by(symbol) %>%
  summarize(n=n()) %>%
  arrange(-n) %>%
  mutate(AMP_tbl = NA)

## HugeAMP gene-associations query ####
return_AMP_genes <- function(gene, index='gene-associations') {
  url <- paste0("https://bioindex.hugeamp.org/api/bio/query/",index,"?q=",gene,"&fmt=row")
  req <- curl_fetch_memory(url)
  nested_list <- rjson::fromJSON(jsonlite::prettify(rawToChar(req$content)))
  AMP_gene_assoc <- map_dfr(nested_list$data, as_tibble)
  return(AMP_gene_assoc)
}
genes2query <- c(topGenes$symbol[1:25],
                 (Genes%>%arrange(minGwasP)%>%filter(term=="PXS_T2D"))$symbol
) %>% unique()
for (i in 1:length(genes2query)) {
  gene <- genes2query[i]
  
  if (exists("AMP_genes")) {
    if (gene %in% AMP_genes$gene) {next}
  }
  
  AMP_tbl <- return_AMP_genes(gene)
  if (nrow(AMP_tbl) == 0) {
    print(paste0("Found no results for ", gene))
    next
  }
  if ("T2D" %in% AMP_tbl$phenotype) {
    T2D_p <- (AMP_tbl %>% filter(phenotype == "T2D"))$pValue
  } else {T2D_p <- as.numeric(NA)}
  
  if (!exists("AMP_genes")) {AMP_genes <- AMP_tbl
  } else {AMP_genes <- AMP_genes %>% add_row(AMP_tbl)}
  
  print(paste0(i," :: p = ",formatC(T2D_p, digits=2, format="E")," :: Queried AMP for ", gene))
  topGenes <- topGenes %>%
    mutate(AMP_tbl = ifelse(symbol == gene, list(.GlobalEnv$AMP_tbl), AMP_tbl))
}

# save to system
#loc_out <- "scratch/AMP_genes.txt"
#fwrite(AMP_genes, loc_out, sep="\t")

# gets AMP T2D gene associations dataset
AMP_T2D_genes <- as_tibble(fread("input_data/AMP_T2D_gene_table.csv"))


AMP_phenos2check <- c("T2D","BMI","TG","HDL","SBP","HBA1C","2hrG","FG","FI")

AMP_summary <- AMP_genes %>%
  select(gene, pValue, phenotype, chromosome, start, end) %>%
  filter(phenotype %in% AMP_phenos2check) %>%
  pivot_wider(names_prefix = "p_",
              names_from = phenotype,
              values_from = pValue) %>%
  select(-p_T2D) %>% # imported manually later
  left_join(
    Genes %>% filter(term=="PXS_T2D") %>% select(symbol, p_PXS_T2D = minGwasP),
    by=c("gene"="symbol") 
  ) %>% full_join(
    AMP_T2D_genes %>% select(gene, p_T2D = minP), by="gene"
  ) %>%
  mutate(p_T2D = ifelse(is.na(p_T2D),1,p_T2D),
         p_PXS_T2D = ifelse(is.na(p_PXS_T2D),1,p_PXS_T2D),
         sig_T2D = p_T2D < 2.5E-6,
         sig_PXS_T2D = p_PXS_T2D < 2.5E-6)

table(AMP_summary %>% select(starts_with("sig_")))

ggplot(AMP_summary, aes(x = log10(-log10(p_T2D)),
                        y = log10(-log10(p_PXS_T2D))
)) +
  geom_point() +
  geom_hline(yintercept = log10(-log10(2.5E-6)), color="red") +
  geom_vline(xintercept = log10(-log10(2.5E-6)), color="red")