library(tidyverse)
library(data.table)
source('paths.R')

dir_FUMA <- paste0(dir_scratch,"FUMA_results/")
# loads job ID table to match jobs to fields
fields <- as_tibble(fread(paste0(dir_scratch,"fields_tbl.txt")))
jobID_tbl <- as_tibble(fread(paste0(dir_FUMA,"jobid_name.csv"))) %>%
  mutate(term = substring(JobName,5)) %>%
  left_join(fields %>% select(term, traitname), by="term") %>%
  select(-V4)
jobID_tbl[jobID_tbl$JobName=="PXS NC",c("term","traitname")] <- list("PXS_T2D","PXS for Type 2 Diabetes")

# looks at GTEx general data
for (i in 1:(nrow(jobID_tbl)-1)) {
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
ggplot(GTEx_general, aes(x=VARIABLE, y=traitname)) +
  geom_tile(aes(fill=-log10(P))) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=30)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=50)) +
  scale_fill_gradient2(low="red",mid="white", high="green", midpoint = -log10(0.05)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Tissue Category") + ylab("Trait") +
  labs(title = "Average gene expression per tissue category for each behavior",
       subtitle = "Using MAGMA with GTEx through FUMA. White: p=0.05 (not adjusted)")

# looks at GTEx specific tissue data
for (i in 1:(nrow(jobID_tbl)-1)) {
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
ggplot(GTEx_specific, aes(x=VARIABLE, y=traitname)) +
  geom_tile(aes(fill=-log10(P))) +
  scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=30)) +
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width=50)) +
  scale_fill_gradient2(low="red",mid="white", high="green", midpoint = -log10(0.05)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Tissue") + ylab("Trait") +
  labs(title = "Average gene expression per tissue for each behavior",
       subtitle = "Using MAGMA with GTEx through FUMA. White: p=0.05 (not adjusted)")

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
