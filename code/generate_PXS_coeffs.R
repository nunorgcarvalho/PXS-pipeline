# exploring T2D phenotype
## Libraries and directories ##
library(tidyverse)
library(data.table)
library(callr)
source('code/paths.R')
library(missMDA)

loc_pheno_full1 <- "/n/no_backup2/patel/uk_biobank/main_data_34521/ukb34521.tab"
loc_pheno_full2 <- "/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv"
loc_40PCs <- "/n/groups/patel/nuno/key_data/UKB_40PCs_500k.txt"

# gets list of exposures
loc_expos <- paste0(dir_script,"../input_data/T2D_XWAS_exposures.txt")
expos_tbl <- as_tibble(fread(locexpos))
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
  # filters out non-British self-ID'd ethnicity
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
fIDs_expos_CRFs <- colnames(data_PHESANT)[2:ncol(data_PHESANT)] %>% str_replace_all("#",".")
fIDs_expos <- fIDs_expos_CRFs[sapply(str_split(fIDs_expos_CRFs,"\\."), `[`,1) %in% expos1]
col_expos <- paste0("f",fIDs_expos)
colnames(data_PHESANT) <- c("userId",paste0("f",fIDs_expos_CRFs))

col_covs <- c("sex", "age", "assessment_center", paste0("pc",1:40))

# impute behavior data ####
data_FAMD <- data_PHESANT %>% select(all_of(col_expos))
FAMD_impute_agency <- missMDA::imputeFAMD(data_FAMD, ncp=4, seed=1) # took about 15-20mins with ncp=4
FAMD_imputed <- FAMD_impute_agency$completeObs
save(FAMD_impute_agency, file='scratch/PHESANT_results/T2D_exposures_FAMD_imputed.RData')

## sets upper and lower limits of imputed data ####
for (term in col_expos) {
  PHESANT_min <- min(data_PHESANT[[term]], na.rm=TRUE)
  PHESANT_max <- max(data_PHESANT[[term]], na.rm=TRUE)
  
  n_under <- sum(FAMD_imputed[[term]] < PHESANT_min)
  if (n_under>0) {print(paste(n_under,'observations under min:',term))}
  FAMD_imputed[FAMD_imputed[[term]] < PHESANT_min, term] <- PHESANT_min
  
  n_over <- sum(FAMD_imputed[[term]] > PHESANT_max)
  if (n_over>0) {print(paste(n_over,'observations over max',term))}
  FAMD_imputed[FAMD_imputed[[term]] > PHESANT_max, term] <- PHESANT_max
}

# makes table of all relevant fields
fields <- tibble(term = c(col_covs,colnames(data_PHESANT[-1]))) %>%
  mutate(field = c(col_covs,fIDs_expos_CRFs)) %>%
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

# makes pheno_EC table that will be saved ####
pheno_EC <- T2D_tbl %>%
  select(FID=userId,IID=userId, T2D_onset=PHENO, T2D_onset_days=TIME,
         sex, age, assessment_center, starts_with("pc")) %>%
  arrange(IID) %>% # PHESANT sorts by ID for some reason
  add_column(FAMD_imputed) %>%
  add_column(data_PHESANT %>% select(all_of(fields$term[fields$use_type=='CRF']))) %>%
  mutate(T2D_onset = as.numeric(T2D_onset))

# defines normal and overweight (0 and 1 respectively) groups
overweight_BMI <- 25
IDs_BMI0 <- T2D_tbl$userId[T2D_tbl$f.21001.0.0 <= overweight_BMI]
IDs_BMI1 <- T2D_tbl$userId[T2D_tbl$f.21001.0.0 > overweight_BMI]
pheno_EC <- pheno_EC %>%
  mutate(overweight = ifelse(IID %in% IDs_BMI0, 0,
                      ifelse(IID %in% IDs_BMI1, 1, NA)))

# assigns testing/training group
set.seed(2024)
pheno_EC$sample_group <- sample(c("A","B"), nrow(pheno_EC), replace=TRUE,
                                prob = c(0.8, 0.2))

# saves pheno_EC to scratch
fwrite(pheno_EC, loc_phenoEC, sep="\t", na="NA", quote=FALSE)

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


# reads pheno_EC back into R ####
# pheno_EC <- as_tibble(fread(loc_phenoEC))

# cross validation cox ridge regression ####


training_tbl <- pheno_EC %>% filter(sample_group=='A')
# assigns equal-sized folds to individuals
K <- 10
N <- nrow(training_tbl)
training_tbl$fold <- ( ceiling( (1:N) / (N/K) ) )[ sample(1:N) ]
#lambdas <- c(0,ncol(pheno_EC))
lambdas <- c(1000,10000,100000)
lambda_iterations <- 10

CV_table <- tibble(lambda=0,k=0,C_stat=0, C_stat_se=0)[0,]
mean_lambda_Cs <- c()

for (i in 1:3) {
  lambda <- lambdas[i]
  
  for (k in 1:K) {
    training_folds <- training_tbl %>% filter(fold!=k)
    testing_fold <- training_tbl %>% filter(fold==k)
    
    training_matrix <- as.matrix( training_folds %>% select(all_of(col_covs[-3]), all_of(col_expos)) )
    testing_matrix <- as.matrix( testing_fold %>% select(all_of(col_covs[-3]), all_of(col_expos)) )
    
    # prepping data entry into coxph function
    ridge.formula <- as.formula(
      paste0("survival::Surv(TIME,PHENO) ~ survival::ridge(training_matrix,theta=",2*lambda,")")
      )
    coxph1 <- survival::coxph(ridge.formula, data = training_folds)
    testing_fold$pred <- scale( testing_matrix %*% coef(coxph1) )[,1]
    
    sc1 <- survival::concordance(survival::Surv(TIME, PHENO) ~ pred, data=testing_fold,
                                 reverse=TRUE)
    
    CV_table <- CV_table %>% add_row(k=k,lambda=lambda,
                                     C_stat = sc1$concordance, C_stat_se = sqrt(sc1$var))
    
    print(paste('lambda =',lambda,':: k =',k,':: C = ', sc1$concordance))
  }
  # mean_lambda_C <- mean( CV_table$C_stat[CV_table$lambda == lambda] )
  # mean_lambda_Cs <- c(mean_lambda_Cs, mean_lambda_C)
  # 
  # if (i==1) {next}
  # 
  # top_lambdas <- lambdas[ order(mean_lambda_Cs, decreasing=TRUE)[1:2] ]
  # next_lambda <- mean(top_lambdas)
  # 
  # if (next_lambda %in% lambdas) {
  #   random_top_lambdas <- lambdas[ order(mean_lambda_Cs, decreasing=TRUE)[1:i] ]
  #   top_lambdas <- sample(random_top_lambdas, 2)
  #   next_lambda <- mean(top_lambdas)
  # }
  # 
  # lambdas <- c(lambdas, next_lambda)
  
}

ggplot(CV_table %>% group_by(lambda) %>%
         summarize(C_stat = mean(C_stat)),
       aes(x=log10(lambda), y=C_stat)) +
  geom_point() +
  # geom_errorbar(aes(ymin=C_stat - C_stat_se,
  #                   ymax=C_stat + C_stat_se)) +
  theme_bw()


loc_CV_table <- 'sandbox/scratch/2024-12_work/CV_table_v1.tsv'
#fwrite(CV_table, loc_CV_table,sep='\t')
#CV_table <- as_tibble(fread(loc_CV_table))










# OLD XWAS + PXS CODE ####

# # randomly sorts individuals into A, B, and C groups
# set.seed(2016)
# pheno_EC$sample_group <- sample(c("A","B","C"), nrow(pheno_EC), replace=TRUE,
#                                 prob = c(0.6, 0.2, 0.2))
# IDA <- (pheno_EC %>% filter(sample_group=="A"))$ID
# IDB <- (pheno_EC %>% filter(sample_group=="B"))$ID
# IDC <- (pheno_EC %>% filter(sample_group=="C"))$ID
# 
# 
# # runs XWAS for both sets of exposures
# xwas1 <- xwas(pheno_EC, X = col_expos, cov = col_covs[-3], mod="cox", IDA = IDA, adjust="fdr") # <-- use if including AC columns
# 
# xwas1_c <- as_tibble(xwas1) %>%
#   mutate(Field = col_expos,
#          FieldID = as.numeric(sapply(str_split(fIDs_expos,"\\."), `[`, 1)),
#          Value = sapply(str_split(fIDs_expos,"\\."), `[`, 2)) %>%
#   mutate(Value = ifelse(Value==100,-7,Value)) %>%
#   left_join(ukb_dict %>% select(FieldID,FieldName=Field, Coding)) %>%
#   left_join(ukb_codings, by=c("Coding","Value")) %>%
#   arrange(-fdr)
# 
# path.out <- 'input_data/xwas_coefficients.txt'
# fwrite(xwas1_c, path.out, sep='\t')
# 
# # calculates the PXS
# sig_expos1_tbl <- (xwas1_c %>%
#                     filter(FieldID %in% (xwas1_c %>% 
#                                            group_by(FieldID) %>% 
#                                            summarize(fdr = min(fdr)) %>%
#                                            filter(fdr < 0.05))$FieldID ) )
# sig_expos1 <- sig_expos1_tbl$Field
# sig_expos1_groups <- sig_expos1_tbl$FieldID
# source(loc_PXS_function)
# PXS1 <- PXS(df = pheno_EC, X = sig_expos1, cov = col_covs[-3], mod = "cox",
#             IDA = IDA, IDB = IDB, IDC = c(), seed = 2016, alph=1)
# 
# PXS1_coeffs <- as_tibble(PXS1) %>%
#   mutate(field = sapply(str_split(term,"f"), `[`, 2)) %>%
#   mutate(fieldID = as.numeric(sapply(str_split(field,"\\."), `[`, 1)),
#          value = sapply(str_split(field,"\\."), `[`, 2)) %>%
#   mutate(value = ifelse(value==100,-7,value)) %>%
#   left_join(ukb_dict %>% select(fieldID=FieldID,fieldname=Field, coding=Coding)) %>%
#   left_join(ukb_codings, by=c("coding"="Coding","value"="Value")) %>%
#   mutate(field = ifelse(is.na(field),term,field),
#          disease = "T2D") %>%
#   select(-std.error,-statistic) %>%
#   arrange(p.value)
# 
# loc_out <- paste0(dir_script,"../input_data/PXS_coefficients.txt")
# fwrite(PXS1_coeffs, loc_out, sep="\t")
