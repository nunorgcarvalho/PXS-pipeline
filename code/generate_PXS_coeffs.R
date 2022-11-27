# exploring T2D phenotype
## Libraries and directories ##
library(tidyverse)
library(data.table)
library(PXStools)
library(callr)

dir_script <- "~/jobs/PXS_pipeline/code/"
dir_scratch <- "~/scratch3/PXS_pipeline/"
loc_pheno_full1 <- "/n/groups/patel/uk_biobank/main_data_34521/ukb34521.tab"
loc_pheno_full2 <- "/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv"
loc_40PCs <- "~/scratch3/key_data/UKB_40PCs_500k.txt"
dir_data_showcase <- "~/scratch3/key_data/" # contains 'Data_Dictionary_Showcase.tsv' from UKBB

# gets list of exposures
loc_expos <- paste0(dir_script,"../input_data/T2D_XWAS_exposures.txt")
expos_tbl <- as_tibble(fread(loc_expos))
expos1 <- (expos_tbl %>% filter(Exposure_Class %in% c("agency")))$FieldID
expos2 <- (expos_tbl %>% filter(Exposure_Class %in% c("agency","no_agency","household")))$FieldID

# gets date of first assessment
cols_to_keep1 <- c("f.eid",paste0("f.",c(53,31,34,54,21000,22021,30750,2976, expos2),".0.0"))
pheno1 <- as_tibble(fread(loc_pheno_full1, select = cols_to_keep1))

# gets date of T2D diagnosis
cols_to_keep2 <- c("eid",paste0(c(130708),"-0.0"))
pheno2 <- as_tibble(fread(loc_pheno_full2, select = cols_to_keep2))

# gets the first 40 principal components
PCs40 <- as_tibble(fread(loc_40PCs)) %>% select(-FID)

# combines both tables
pheno <- pheno1 %>%
  left_join(pheno2, by=c("f.eid"="eid")) %>%
  left_join(PCs40, by=c("f.eid"="IID"))

# defines whether an individual got T2D after their first assessment
T2D_tbl <- pheno %>%
  rename(IID = f.eid, assessment1_date = `f.53.0.0`, T2D_date = `130708-0.0`, ethnicity = `f.21000.0.0`,
         HbA1c = `f.30750.0.0`, sex = `f.31.0.0`, age = `f.34.0.0`, assessment_center = `f.54.0.0`, T2D_report_age = `f.2976.0.0`) %>%
  filter(ethnicity == 1001, #`f.22021.0.0` == 0,
         HbA1c < 48) %>%
  filter_at(vars(starts_with("pc"), sex, age, assessment_center), all_vars(!is.na(.))) %>%
  mutate(T2D_date = ifelse(T2D_date < "2037-07-07" & T2D_date > "1903-03-03",T2D_date, NA)) %>%
  mutate(T2D_date = as.Date.numeric(T2D_date, "1970-01-01")) %>%
  mutate(T2D_after_assessment = ifelse(!is.na(T2D_date),T2D_date > assessment1_date, FALSE),
         T2D_before_assessment = ifelse(!is.na(T2D_date),T2D_date <= assessment1_date, FALSE),
         T2D_reported = (T2D_report_age > 0) & !(is.na(T2D_report_age))) %>%
  mutate(T2D_onset = ifelse(T2D_after_assessment,as.numeric(T2D_date - assessment1_date),max(as.numeric(T2D_date - assessment1_date), na.rm=TRUE))) %>%
  filter(!T2D_before_assessment, !T2D_reported) %>%
  rename(userId = IID, PHENO = T2D_after_assessment, TIME = T2D_onset)

# number of cases
sum(T2D_tbl$PHENO)

# makes table formatted for PHESANT to standardize exposures
tbl_out <- T2D_tbl[,c("userId",paste0("f.",expos2,".0.0"))]
colnames(tbl_out) <- c("userId",paste0("x",expos2,"_0_0"))
loc_out <- paste0(dir_scratch,"T2D_exposures_tbl_raw.txt")
write.table(tbl_out, loc_out, sep="\t", quote=FALSE, row.names=FALSE)

# RUN PHESANT THROUGH 'run_PHESANT.sh'
# this R code runs 'Rscript' from within R ratehr than the CMD
# takes a few minutes to run
dir_PHESANT_WAS <- "~/PHESANT/WAS/"
rscript(script="phenomeScan.r",
        cmdargs = c(paste0('--phenofile=',dir_scratch,'T2D_exposures_tbl_raw.txt'),
                    paste0('--variablelistfile=',dir_PHESANT_WAS,'../variable-info/outcome-info.tsv'),
                    paste0('--datacodingfile=',dir_PHESANT_WAS,'../variable-info/data-coding-ordinal-info.txt'),
                    paste0('--resDir=',dir_scratch,'/PHESANT_results/'),
                    '--genetic=FALSE',
                    '--save',
                    '--tab=TRUE'),
        echo = TRUE, wd = dir_PHESANT_WAS)

# load PHESANT results
data_cont <- as_tibble(fread(paste0(dir_scratch,"PHESANT_results/data-cont-all.txt")))
data_catunord <- as_tibble(fread(paste0(dir_scratch,"PHESANT_results/data-catunord-all.txt")))
data_catord <- as_tibble(fread(paste0(dir_scratch,"PHESANT_results/data-catord-all.txt")))
data_binary <- as_tibble(fread(paste0(dir_scratch,"PHESANT_results/data-binary-all.txt")))

data_PHESANT <- data_cont %>% left_join(data_catunord) %>% left_join(data_catord) %>% left_join(data_binary)
PHESANT_fIDs2 <- colnames(data_PHESANT)[2:ncol(data_PHESANT)] %>% str_replace_all("#",".")
PHESANT_fIDs1 <- PHESANT_fIDs2[sapply(str_split(PHESANT_fIDs2,"\\."), `[`,1) %in% expos1]
col_expos1 <- paste0("f",PHESANT_fIDs1)
col_expos2 <- paste0("f",PHESANT_fIDs2)
colnames(data_PHESANT) <- c("userId",col_expos2)

T2D_tbl2 <- T2D_tbl %>%
  select(userId, PHENO, TIME, sex, age, assessment_center, starts_with("pc")) %>%
  left_join(data_PHESANT) %>%
  rename(ID = userId)

# randomly sorts individuals into A, B, and C groups
set.seed(2016)
T2D_tbl2$sample_group <- sample(c("A","B","C"), nrow(T2D_tbl2), replace=TRUE)
#T2D_tbl2 <- T2D_tbl2_full %>% filter(row_number() %in% sample(1:nrow(T2D_tbl2), 20000))
T2D_tbl2 <- T2D_tbl2_full
IDA <- (T2D_tbl2 %>% filter(sample_group=="A"))$ID
IDB <- (T2D_tbl2 %>% filter(sample_group=="B"))$ID
IDC <- (T2D_tbl2 %>% filter(sample_group=="C"))$ID

#subset <- T2D_tbl2 %>% filter(row_number() %in% sample(1:nrow(T2D_tbl), 20000))

# runs XWAS for both sets of exposures
col_covs <- c("sex", "age", "assessment_center", paste0("pc",1:40))
xwas1 <- xwas(T2D_tbl2, X = col_expos1, cov = col_covs, mod="cox", IDA = IDA, adjust="fdr")
xwas2 <- xwas(T2D_tbl2, X = col_expos2, cov = col_covs, mod="cox", IDA = IDA, adjust="fdr")


ukb_dict <- as_tibble(fread(paste0(dir_data_showcase,"Data_Dictionary_Showcase.tsv")))
ukb_codings <- as_tibble(fread(paste0(dir_data_showcase,"Codings.tsv"),quote=""))
xwas1_c <- as_tibble(xwas1) %>%
  mutate(Field = col_expos1,
         FieldID = as.numeric(sapply(str_split(PHESANT_fIDs1,"\\."), `[`, 1)),
         Value = sapply(str_split(PHESANT_fIDs1,"\\."), `[`, 2)) %>%
  mutate(Value = ifelse(Value==100,-7,Value)) %>%
  left_join(ukb_dict %>% select(FieldID,FieldName=Field, Coding)) %>%
  left_join(ukb_codings, by=c("Coding","Value")) %>%
  arrange(-fdr)
xwas2_c <- as_tibble(xwas2) %>%
  mutate(Field = col_expos2,
         FieldID = as.numeric(sapply(str_split(PHESANT_fIDs2,"\\."), `[`, 1)),
         Value = sapply(str_split(PHESANT_fIDs2,"\\."), `[`, 2)) %>%
  mutate(Value = ifelse(Value==100,-7,Value)) %>%
  left_join(ukb_dict %>% select(FieldID,FieldName=Field, Coding)) %>%
  left_join(ukb_codings, by=c("Coding","Value")) %>%
  arrange(-fdr)

# calculates the PXS
sig_expos1 <- (xwas1_c %>% filter(fdr < 0.05))$Field
sig_expos2 <- (xwas2_c %>% filter(fdr < 0.05))$Field
library(loggr)
library(glmnet)
#source(paste0(dir_script,"PXS_modified.R"))
source(paste0(dir_script,"../../PXStools/R/PXS.R"))
PXS1 <- PXS(df = T2D_tbl2, X = sig_expos1, cov = col_covs, mod = "cox",
            IDA = c(IDA,IDC), IDB = IDB, IDC = c(), seed = 2016, alph=1)
PXS1 <- PXS_modified(df = T2D_tbl2, X = sig_expos1, cov = col_covs, mod = "cox",
                     IDA = c(IDA, IDC), IDB = IDB, IDC = c(), seed = 2016, alph=1)

PXS1_coeffs <- as_tibble(PXS1) %>%
  mutate(Field = sapply(str_split(term,"f"), `[`, 2)) %>%
  mutate(FieldID = as.numeric(sapply(str_split(Field,"\\."), `[`, 1)),
         Value = sapply(str_split(Field,"\\."), `[`, 2)) %>%
  mutate(Value = ifelse(Value==100,-7,Value)) %>%
  left_join(ukb_dict %>% select(FieldID,FieldName=Field, Coding)) %>%
  left_join(ukb_codings, by=c("Coding","Value")) %>%
  mutate(Field = ifelse(is.na(Field),term,Field)) %>%
  select(-term, -std.error,-statistic) %>%
  arrange(p.value)

PXS2 <- PXS(df = T2D_tbl2, X = sig_expos2, cov = col_covs, mod = "cox",
            IDA = c(IDA, IDC), IDB = IDB, IDC = c(), seed = 2016, alph=1)

PXS2_coeffs <- as_tibble(PXS2) %>%
  mutate(Field = sapply(str_split(term,"f"), `[`, 2)) %>%
  mutate(FieldID = as.numeric(sapply(str_split(Field,"\\."), `[`, 1)),
         Value = sapply(str_split(Field,"\\."), `[`, 2)) %>%
  mutate(Value = ifelse(Value==100,-7,Value)) %>%
  left_join(ukb_dict %>% select(FieldID,FieldName=Field, Coding)) %>%
  left_join(ukb_codings, by=c("Coding","Value")) %>%
  mutate(Field = ifelse(is.na(Field),term,Field)) %>%
  select(-term, -std.error,-statistic) %>%
  arrange(p.value)

loc_out <- paste0(dir_script,"../input_data/T2D_PXS_coefficients.txt")
write.table(PXS1_coeffs, loc_out, sep="\t", row.names=FALSE, quote=FALSE)
