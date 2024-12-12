library(tidyverse)
library(data.table)
library(callr)
source('code/paths.R')
source('../ePC/code/helper_functions.R')

# Load ####
raw_phenofiles <- c("/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv",
                    "/n/no_backup2/patel/uk_biobank/main_data_34521/ukb34521.csv",
                    "/n/no_backup2/patel/uk_biobank/main_data_9512/ukb9512.csv")
loc_40PCs <- "/n/groups/patel/nuno/key_data/UKB_40PCs_500k.txt"
loc_expos <- paste0(dir_script,"../input_data/T2D_XWAS_exposures.txt")
expos_tbl <- as_tibble(fread(loc_expos))
expos1 <- (expos_tbl %>% filter(Exposure_Class %in% c("agency")))$FieldID


cols_to_keep <- c(53,31,34,54,21000,22021,40000, expos1) %>% unique()
pheno_raw <- extract_phenos(cols_to_keep, raw_phenofiles)
colnames(pheno_raw)[-1] <- str_replace(colnames(pheno_raw)[-1],'-0.0','')

tbl_out <- pheno_raw %>% select(userId='eid', as.character(expos1))
colnames(tbl_out)[-1] <- paste0('x',colnames(tbl_out)[-1],'_0_0')

# PHESANT ####
dir.create(paste0('scratch/mortality/PHESANT_results/'))
loc_out <- paste0('scratch/mortality/PHESANT_results/death_exposures_tbl_raw.txt')
fwrite(tbl_out, loc_out, sep="\t")

# this R code runs 'Rscript' from within R rather than the CMD
# takes a few minutes to run
dir_PHESANT_WAS <- "~/PHESANT/WAS/"
rscript(script="phenomeScan.r",
        cmdargs = c(paste0('--phenofile=',dir_scratch,'mortality/PHESANT_results/death_exposures_tbl_raw.txt'),
                    paste0('--variablelistfile=',dir_PHESANT_WAS,'../variable-info/outcome-info.tsv'),
                    paste0('--datacodingfile=',dir_PHESANT_WAS,'../variable-info/data-coding-ordinal-info.txt'),
                    paste0('--resDir=',dir_scratch,'mortality/PHESANT_results/'),
                    '--genetic=FALSE',
                    '--save',
                    '--tab=TRUE'),
        echo = TRUE, wd = dir_PHESANT_WAS)

# load PHESANT results
data_cont <- as_tibble(fread(paste0(dir_scratch,"mortality/PHESANT_results/data-cont-all.txt")))
data_catunord <- as_tibble(fread(paste0(dir_scratch,"mortality/PHESANT_results/data-catunord-all.txt")))
data_catord <- as_tibble(fread(paste0(dir_scratch,"mortality/PHESANT_results/data-catord-all.txt")))
data_binary <- as_tibble(fread(paste0(dir_scratch,"mortality/PHESANT_results/data-binary-all.txt")))


# fixes catunord data to be binary multi-column table instead
data_catunord2 <- tibble(userID = data_catunord$userID)
for (col in colnames(data_catunord)[-1]) {
  print(col)
  values <- data_catunord[,col] %>% drop_na() %>% unlist() %>% unique()
  for (value in values) {
    col_value <- paste0(col,"#",value)
    data_catunord2[,col_value] <- as.numeric(data_catunord[,col] == value)
  }
}
# changes binary columns coded as 1,2 to be 0,1
data_binary2 <- tibble(userID = data_binary$userID)
for (col in colnames(data_binary)[-1]) {
  print(col)
  values <- data_binary[,col] %>% drop_na() %>% unlist() %>% unique()
  if (sort(values) == c(1,2)) {data_binary2[,col] <- data_binary[,col] - 1
  } else {data_binary2[,col] <- data_binary[,col]}
}

data_PHESANT <- data_cont %>% left_join(data_catunord2) %>% left_join(data_catord) %>% left_join(data_binary2)
PHESANT_fIDs <- colnames(data_PHESANT)[2:ncol(data_PHESANT)] %>% str_replace_all("#",".")
col.expos <- paste0("f",PHESANT_fIDs)
colnames(data_PHESANT) <- c("userId",paste0("f",PHESANT_fIDs))

# make big table ####
PCs40 <- as_tibble(fread(loc_40PCs)) %>% select(-FID)

pheno <- pheno_raw %>%
  select(FID = eid, IID=eid, sex=`31`, yob=`34`, assessment_center=`54`,
         ethnicity=`21000`,ass1.date=`53`, death.date=`40000`) %>%
  left_join(data_PHESANT, by=c('IID'='userId')) %>%
  left_join(PCs40, by='IID') %>%
  mutate(died = !is.na(death.date),
         death_ass1.diff = ifelse(died,
                                  as.numeric(death.date) - as.numeric(ass1.date),
                                  max(as.numeric(death.date) - as.numeric(ass1.date),
                                      na.rm=TRUE))
  )

loc.out <- 'scratch/mortality/pheno_EC.death.txt'
fwrite(pheno, loc.out,sep='\t')

# XWAS ####
col.covs <- c("sex", "yob", "assessment_center", paste0("pc",1:40))

set.seed(2016)
pheno_PXS <- pheno %>%
  rename(PHENO = died, TIME = death_ass1.diff, ID = IID) %>%
  mutate(sample_group = sample(c("A","B","C"), nrow(pheno), replace=TRUE)) %>%
  filter(ethnicity == '1001')

IDA <- (pheno_PXS %>% filter(sample_group=="A"))$ID
IDB <- (pheno_PXS %>% filter(sample_group=="B"))$ID
IDC <- (pheno_PXS %>% filter(sample_group=="C"))$ID

library(PXStools)
xwas1 <- xwas(pheno_PXS, X = col.expos, cov = col.covs[-3], mod="cox", IDA = c(IDA,IDC), adjust="fdr")

# reads UKB field data
fields_tbl <- as_tibble(fread('scratch/fields_tbl.txt'))
ukb_dict <- as_tibble(fread(paste0(dir_data_showcase,"Data_Dictionary_Showcase.tsv")))
ukb_codings <- as_tibble(fread(paste0(dir_data_showcase,"Codings.tsv"),quote=""))


xwas1_c <- as_tibble(xwas1) %>%
  mutate(Field = col.expos,
         FieldID = as.numeric(sapply(str_split(PHESANT_fIDs,"\\."), `[`, 1)),
         Value = sapply(str_split(PHESANT_fIDs,"\\."), `[`, 2)) %>%
  mutate(Value = ifelse(Value==100,-7,Value)) %>%
  left_join(ukb_dict %>% select(FieldID,FieldName=Field, Coding)) %>%
  left_join(ukb_codings, by=c("Coding","Value")) %>%
  arrange(-fdr)

path.out <- 'scratch/mortality/xwas_coefficients.txt'
fwrite(xwas1_c, path.out, sep='\t')

sig_expos1.tbl <- (xwas1_c %>%
                     filter(FieldID %in% (xwas1_c %>% 
                                            group_by(FieldID) %>% 
                                            summarize(fdr = min(fdr)) %>%
                                            filter(fdr < 0.05))$FieldID ) )
sig_expos1 <- sig_expos1.tbl$Field
sig_expos1.groups <- sig_expos1.tbl$FieldID


# PXS ####
loc_PXS_function <- paste0(dir_script,"../../PXStools/R/PXS.R")
source(loc_PXS_function)

PXS1 <- PXS(df = pheno_PXS, X = sig_expos1, cov = col.covs[-3], mod = "cox",
            IDA = c(IDA,IDC), IDB = IDB, IDC = c(), seed = 2016, alph=1)

PXS1_coeffs <- as_tibble(PXS1) %>%
  mutate(field = sapply(str_split(term,"f"), `[`, 2)) %>%
  mutate(fieldID = as.numeric(sapply(str_split(field,"\\."), `[`, 1)),
         value = sapply(str_split(field,"\\."), `[`, 2)) %>%
  mutate(value = ifelse(value==100,-7,value)) %>%
  left_join(ukb_dict %>% select(fieldID=FieldID,fieldname=Field, coding=Coding)) %>%
  left_join(ukb_codings, by=c("coding"="Coding","value"="Value")) %>%
  mutate(field = ifelse(is.na(field),term,field),
         disease = "death") %>%
  select(-std.error,-statistic) %>%
  arrange(p.value)

loc_out <- paste0(dir_script,"../scratch/mortality/PXS_coefficients.txt")
fwrite(PXS1_coeffs, loc_out, sep="\t")
