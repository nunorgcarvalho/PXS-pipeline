# exploring T2D phenotype
## Libraries and directories ##
library(tidyverse)
library(data.table)
library(callr)
source('code/paths.R')
library(PXStools)

loc_pheno_full1 <- "/n/no_backup2/patel/uk_biobank/main_data_34521/ukb34521.tab"
loc_pheno_full2 <- "/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv"
loc_40PCs <- "/n/groups/patel/nuno/key_data/UKB_40PCs_500k.txt"
# location to modified PXS function as present here:
# https://github.com/nunorgcarvalho/PXStools/blob/nuno_edits/R/PXS.R
loc_PXS_function <- paste0(dir_script,"../../PXStools/R/PXS.R")

# gets list of exposures
loc_expos <- paste0(dir_script,"../input_data/T2D_XWAS_exposures.txt")
expos_tbl <- as_tibble(fread(loc_expos))
expos1 <- (expos_tbl %>% filter(Exposure_Class %in% c("agency")))$FieldID

# loads list of CRFs
loc_CFRs_tbl <- paste0(dir_script,"../input_data/CRFs_table.txt")
CRFs_tbl <- as_tibble(fread(loc_CFRs_tbl))
CRFs_tbl <- CRFs_tbl %>%
  mutate(fieldID = sapply(str_split(CRFs_tbl$field,"\\."), `[`,2))

# gets date of first assessment
cols_to_keep1 <- c("f.eid",paste0("f.",c(53,31,34,54,21000,22021,30750,2976, expos1),".0.0"), CRFs_tbl$field_raw) %>% unique()
pheno1 <- as_tibble(fread(loc_pheno_full1, select = cols_to_keep1))

# gets date of T2D diagnosis
cols_to_keep2 <- c("eid",paste0(c(130708),"-0.0"), CRFs_tbl$field_raw)
pheno2 <- as_tibble(fread(loc_pheno_full2, select = cols_to_keep2))

# gets the first 40 principal components
PCs40 <- as_tibble(fread(loc_40PCs)) %>% select(-FID)

# reads UKB field data
ukb_dict <- as_tibble(fread(paste0(dir_data_showcase,"Data_Dictionary_Showcase.tsv")))
ukb_codings <- as_tibble(fread(paste0(dir_data_showcase,"Codings.tsv"),quote=""))


# combines both tables
pheno <- pheno1 %>%
  left_join(pheno2, by=c("f.eid"="eid")) %>%
  left_join(PCs40, by=c("f.eid"="IID"))

colnames(pheno)[colnames(pheno) %in% CRFs_tbl$field_raw] <- sapply(colnames(pheno)[colnames(pheno) %in% CRFs_tbl$field_raw],
                                                                   function(x) (CRFs_tbl%>%filter(field_raw==x))$field)

# defines whether an individual got T2D and when
T2D_tbl <- pheno %>%
  # renames variables
  rename(IID=f.eid,assessment1_date=`f.53.0.0`,T2D_date=`130708-0.0`,
         ethnicity=`f.21000.0.0`,sex=`f.31.0.0`,age=`f.34.0.0`,
         assessment_center=`f.54.0.0`,TxD_report_age = `f.2976.0.0`) %>%
  # filters out non-British ethnicity
  filter(ethnicity == 1001) %>%
  # removes individuals with any missing covariate data
  filter_at(vars(starts_with("pc"), sex, age, assessment_center), all_vars(!is.na(.))) %>%
  # defines T2D_date as valid date of E11 ICD10 reporting
  mutate(T2D_date = ifelse(T2D_date < "2037-07-07" & T2D_date > "1903-03-03",T2D_date, NA)) %>%
  mutate(T2D_date = as.Date.numeric(T2D_date, "1970-01-01")) %>%
  # defines whether an individual had E11 reported after first assessment
  mutate(T2D_after_assessment = ifelse(!is.na(T2D_date),T2D_date > assessment1_date, FALSE),
         T2D_before_assessment = ifelse(!is.na(T2D_date),T2D_date <= assessment1_date, FALSE),
         TxD_reported = (TxD_report_age > 0) & !(is.na(TxD_report_age))) %>%
  # defines days between E11 report and first assessment
  mutate(T2D_onset_days = ifelse(T2D_after_assessment,
                                 as.numeric(T2D_date - assessment1_date),
                                 max(as.numeric(T2D_date - assessment1_date),
                                     na.rm=TRUE)),
         # defines whether individual ever had E11 report
         T2D_all = (T2D_after_assessment | T2D_before_assessment),
         # defines whether individual ever had E11 report or diabetes diagnosis (by doctor or by HbA1c)
         TxD_all = (T2D_after_assessment | T2D_before_assessment |
                    TxD_reported | (`f.30750.0.0` >= 48 & !is.na(`f.30750.0.0`)))) %>%
  # tweaked T2D after first assessment definition that NA's pre-assessment diabetics
  mutate(T2D_onset = ifelse((TxD_all & !T2D_after_assessment) | TxD_reported | `f.30750.0.0` >= 48,
                            NA, T2D_after_assessment))

# keeps table with all the T2D definitions with expanded cohort size. Later sized down
T2D_definitions <- T2D_tbl

# removes pre-assessment diabetics from main table
T2D_tbl <- T2D_tbl %>%
  filter(!is.na(T2D_onset)) %>%
  rename(userId = IID, PHENO = T2D_onset, TIME = T2D_onset_days)

# number of cases
sum(T2D_tbl$PHENO)

# makes table formatted for PHESANT to standardize exposures
tbl_out <- T2D_tbl[,c("userId",paste0("f.",c(expos1,CRFs_tbl$fieldID),".0.0"))]
colnames(tbl_out) <- c("userId",paste0("x",c(expos1,CRFs_tbl$fieldID),"_0_0"))
dir.create(paste0(dir_scratch,'PHESANT_results/'))
loc_out <- paste0(dir_scratch,"PHESANT_results/T2D_exposures_tbl_raw.txt")
fwrite(tbl_out, loc_out, sep="\t")

# this R code runs 'Rscript' from within R rather than the CMD
# takes a few minutes to run
dir_PHESANT_WAS <- "~/PHESANT/WAS/"
rscript(script="phenomeScan.r",
        cmdargs = c(paste0('--phenofile=',dir_scratch,'PHESANT_results/T2D_exposures_tbl_raw.txt'),
                    paste0('--variablelistfile=',dir_PHESANT_WAS,'../variable-info/outcome-info.tsv'),
                    paste0('--datacodingfile=',dir_PHESANT_WAS,'../variable-info/data-coding-ordinal-info.txt'),
                    paste0('--resDir=',dir_scratch,'PHESANT_results/'),
                    '--genetic=FALSE',
                    '--save',
                    '--tab=TRUE'),
        echo = TRUE, wd = dir_PHESANT_WAS)

# load PHESANT results
data_cont <- as_tibble(fread(paste0(dir_scratch,"PHESANT_results/data-cont-all.txt")))
data_catunord <- as_tibble(fread(paste0(dir_scratch,"PHESANT_results/data-catunord-all.txt")))
data_catord <- as_tibble(fread(paste0(dir_scratch,"PHESANT_results/data-catord-all.txt")))
data_binary <- as_tibble(fread(paste0(dir_scratch,"PHESANT_results/data-binary-all.txt")))

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
PHESANT_fIDs2 <- colnames(data_PHESANT)[2:ncol(data_PHESANT)] %>% str_replace_all("#",".")
PHESANT_fIDs1 <- PHESANT_fIDs2[sapply(str_split(PHESANT_fIDs2,"\\."), `[`,1) %in% expos1]
col_expos1 <- paste0("f",PHESANT_fIDs1)
colnames(data_PHESANT) <- c("userId",paste0("f",PHESANT_fIDs2))

T2D_tbl2 <- T2D_tbl %>%
  select(userId, PHENO, TIME, sex, age, assessment_center, starts_with("pc")) %>%
  left_join(data_PHESANT) %>%
  rename(ID = userId)

####
# Makes a column for each assessment_center
# Currently skipping, does not significantly impact results

# ACs <- T2D_tbl2$assessment_center %>% unique()
# for (AC in ACs) {
#   col_AC <- paste0("assessment_center.",AC)
#   T2D_tbl2[,col_AC] <- as.numeric(T2D_tbl2$assessment_center == AC)
# }
####

col_covs <- c("sex", "age", "assessment_center", paste0("pc",1:40))
#col_covs <- c("sex", "age", paste0("assessment_center.",ACs), paste0("pc",1:40)) # <-- use if including AC columns
#fields <- tibble(term = colnames(T2D_tbl2)[-c(1:3,6)]) %>% # <-- use if including AC columns
fields <- tibble(term = colnames(T2D_tbl2)[-c(1:3)]) %>%
  mutate(field = c(col_covs,PHESANT_fIDs2)) %>%
  mutate(use_type = ifelse(field %in% col_covs, "covar",
                           ifelse(field %in% CRFs_tbl$fieldID, "CRF","exposure")),
         data_type = ifelse(str_replace(field,"\\.","#") %in% colnames(data_binary2), "binary",
                     ifelse(str_replace(field,"\\.","#") %in% colnames(data_catord), "catord",
                     ifelse(str_replace(field,"\\.","#") %in% colnames(data_catunord2), "catunord",
                     ifelse(str_replace(field,"\\.","#") %in% colnames(data_cont), "continuous",as.character(NA))))),
         fieldID = sapply(str_split(field,"\\."), `[`, 1),
         value = sapply(str_split(field,"\\."), `[`, 2)) %>%
  mutate(fieldID = as.numeric(ifelse(use_type=="covar",NA,fieldID)),
         value = ifelse(value==100,-7,value)) %>%
  left_join(ukb_dict %>% select(fieldID=FieldID,fieldname=Field, coding=Coding)) %>%
  left_join(ukb_codings, by=c("coding"="Coding","value"="Value")) %>%
  mutate(traitname = ifelse(is.na(Meaning),fieldname,
                            paste0(fieldname,": ",Meaning)))

loc_out <- paste0(dir_scratch,"fields_tbl.txt")
fwrite(fields, loc_out, sep="\t")

tbl_out <- as_tibble(cbind(FID = T2D_tbl2$ID,IID = T2D_tbl2$ID,
                           T2D_tbl2[c(4:ncol(T2D_tbl2),2,3)])) %>%
  rename(T2D_onset = PHENO, T2D_onset_days = TIME) %>%
  mutate(T2D_onset = as.numeric(T2D_onset))
#loc_out <- paste0(dir_scratch,"pheno_EC.txt")
# loc_phenoEC comes from paths.R
fwrite(tbl_out, loc_phenoEC, sep="\t", na="NA", quote=FALSE)

# saves the phenotype data with the T2D definitions
CRFs_tbl <- CRFs_tbl %>% mutate(term = paste0("f",fieldID))
T2D_definitions_out <- tibble(.rows = nrow(T2D_definitions))
# I am ashamed of this code, there has got to be a better way to do this
for (col in colnames(T2D_definitions)) {
  if (col=="IID") {
    T2D_definitions_out[,c("FID","IID")] <- cbind(T2D_definitions$IID, T2D_definitions$IID)
  } else if (col %in% col_covs) {
    T2D_definitions_out[,col] <- T2D_definitions[,col]
  } else if (col %in% CRFs_tbl$field) {
    col2 <- CRFs_tbl$term[CRFs_tbl$field==col]
    T2D_definitions_out[,col2] <- T2D_definitions[,col]
  } else if (grepl("T2D",col) | grepl("TxD",col)) {
    T2D_definitions_out[,col] <- T2D_definitions[,col]
  }
}
loc_out <- paste0(dir_scratch, "phenoEC_fullT2D.txt")
fwrite(T2D_definitions_out, loc_out, sep="\t", na="NA", quote=FALSE, logical01=TRUE)

# some of the following code is commented out because we are opting to use the
# exposures in group 1, not group 2

# randomly sorts individuals into A, B, and C groups
set.seed(2016)
T2D_tbl2$sample_group <- sample(c("A","B","C"), nrow(T2D_tbl2), replace=TRUE,
                                prob = c(0.6, 0.2, 0.2))
# T2D_tbl2 <- T2D_tbl2_full %>% filter(row_number() %in% sample(1:nrow(T2D_tbl2), 20000))
# T2D_tbl2 <- T2D_tbl2_full
IDA <- (T2D_tbl2 %>% filter(sample_group=="A"))$ID
IDB <- (T2D_tbl2 %>% filter(sample_group=="B"))$ID
IDC <- (T2D_tbl2 %>% filter(sample_group=="C"))$ID


# runs XWAS for both sets of exposures
xwas1 <- xwas(T2D_tbl2, X = col_expos1, cov = col_covs[-3], mod="cox", IDA = IDA, adjust="fdr") # <-- use if including AC columns
#xwas1 <- xwas(T2D_tbl2, X = col_expos1, cov = col_covs, mod="cox", IDA = c(IDA,IDC), adjust="fdr")

xwas1_c <- as_tibble(xwas1) %>%
  mutate(Field = col_expos1,
         FieldID = as.numeric(sapply(str_split(PHESANT_fIDs1,"\\."), `[`, 1)),
         Value = sapply(str_split(PHESANT_fIDs1,"\\."), `[`, 2)) %>%
  mutate(Value = ifelse(Value==100,-7,Value)) %>%
  left_join(ukb_dict %>% select(FieldID,FieldName=Field, Coding)) %>%
  left_join(ukb_codings, by=c("Coding","Value")) %>%
  arrange(-fdr)

path.out <- 'input_data/xwas_coefficients.txt'
fwrite(xwas1_c, path.out, sep='\t')

# calculates the PXS
#sig_expos1 <- (xwas1_c %>% filter(fdr < 0.05))$Field
sig_expos1_tbl <- (xwas1_c %>%
                    filter(FieldID %in% (xwas1_c %>% 
                                           group_by(FieldID) %>% 
                                           summarize(fdr = min(fdr)) %>%
                                           filter(fdr < 0.05))$FieldID ) )
sig_expos1 <- sig_expos1_tbl$Field
sig_expos1_groups <- sig_expos1_tbl$FieldID
source(loc_PXS_function)
PXS1 <- PXS(df = T2D_tbl2, X = sig_expos1, cov = col_covs[-3], mod = "cox",
#PXS1 <- PXS(df = T2D_tbl2, X = sig_expos1, cov = col_covs, mod = "cox",
            IDA = IDA, IDB = IDB, IDC = c(), seed = 2016, alph=1)

PXS1_coeffs <- as_tibble(PXS1) %>%
  mutate(field = sapply(str_split(term,"f"), `[`, 2)) %>%
  mutate(fieldID = as.numeric(sapply(str_split(field,"\\."), `[`, 1)),
         value = sapply(str_split(field,"\\."), `[`, 2)) %>%
  mutate(value = ifelse(value==100,-7,value)) %>%
  left_join(ukb_dict %>% select(fieldID=FieldID,fieldname=Field, coding=Coding)) %>%
  left_join(ukb_codings, by=c("coding"="Coding","value"="Value")) %>%
  mutate(field = ifelse(is.na(field),term,field),
         disease = "T2D") %>%
  select(-std.error,-statistic) %>%
  arrange(p.value)

loc_out <- paste0(dir_script,"../input_data/PXS_coefficients.txt")
fwrite(PXS1_coeffs, loc_out, sep="\t")
