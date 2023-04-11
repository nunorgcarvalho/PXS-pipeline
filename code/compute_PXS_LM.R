disease <- commandArgs(trailingOnly = TRUE)[1]
sink(paste0(disease,"_compute_PXS_LM.Rout"))

## Libraries and directories ##
library(tidyverse)
library(data.table)

dir_script <- "~/jobs/PXS_pipeline/code/"
dir_data_showcase <- "~/scratch3/key_data/" # contains 'Data_Dictionary_Showcase.tsv' and 'Codings.tsv' from UKBB

## Code ###
dir_out <- "./"

#col_coeff <- paste0("estimate_",disease)
col_coeff <- "estimate"

# Loads the phenotype file created in initial_setup.R
loc_ukbpheno <- "../pheno_EC.txt"
pheno <- as_tibble(fread(loc_ukbpheno))
loc_fields <- "../fields_tbl.txt"
# fields <- as_tibble(fread("~/scratch3/PXS_pipeline/fields_tbl.txt"))
fields <- as_tibble(fread(loc_fields)) #%>%
#  filter((!!as.name(disease) != 0) | (use_type=="covar"))
col_covs <- (fields %>% filter(use_type=="covar"))$term
#fields <- fields %>% filter(use_type=="exposure")

# loads the file containing XWAS coefficients and filters to just to the disease
loc_coeffs <- paste0(dir_script,"../input_data/PXS_coefficients.txt")
coeffs <- as_tibble(fread(loc_coeffs)) %>%
  select(term, estimate, disease) #%>%
#  filter(!is.na(!!as.name(col_coeff))) %>% select(V1,X,term,ends_with(disease))

# Loads the 'Codings.tsv' file and 'Data_Dictionary_Showcase.tsv' files from
# https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide
# ukb_codings <- as_tibble(fread(paste0(dir_data_showcase,"Codings.tsv"),quote=""))
# ukb_dict <- as_tibble(fread(paste0(dir_data_showcase,"Data_Dictionary_Showcase.tsv")))


# determines which exposures are categorical variables
# categorical_exposures <- (fields %>% filter(use_type=="exposure",
#                                             data_type=="categorical"))$field

# splits the field column in the coeffs file based on the actual field and the
# answer to the prompt for categorical variables
# splits <- as_tibble(str_split_fixed(coeffs$term,regex("\\.0\\.0"),n=2))
# coeffs$field <- paste0(splits$V1,".0.0")
# coeffs$coding_text <- splits$V2
# coeffs$coding <- ""
# coeffs <- coeffs %>% select(X,field, coding, coding_text, ends_with(disease))

# converts the text answer of categorical variables to the numeric answer used
# in the actual data
# for (i in 1:nrow(coeffs)) {
#   the_field <- coeffs$field[i]
#   the_coding_text <- coeffs$coding_text[i]
#   field_id <- coeffs$X[i]
#   
#   if ((nchar(the_coding_text) > 1) & (the_field %in% categorical_exposures))  {
#     print(paste(the_field,the_coding_text))
#     
#     coding_num <- (ukb_dict %>% filter(FieldID==field_id))$Coding
#     coding <- (ukb_codings %>% filter(Coding == coding_num,
#                                       Meaning == the_coding_text)
#     )$Value
#     coeffs$coding[i] <- as.character(coding)
#   }
# }

# # establishes which fields will be used in analysis
# fields <- fields %>% mutate(use_in_PXS = FALSE) #, code = as.character(code))
# for (i in 1:nrow(fields)) {
#   field <- fields$field[i]
#   if ( (field %in% colnames(pheno)) && (fields$use_type[i]=="exposure") ) {
#     fields$use_in_PXS[i] <- TRUE
#   }
# }
# col_use_in_PXS <- paste0("use_in_PXS_",disease)
# colnames(fields)[which(colnames(fields) == "use_in_PXS")] <- col_use_in_PXS
# 
# filters list to just exposures for easy dot-product computation
# exposures_tbl <- fields %>% full_join(coeffs, by="term") %>%
#   filter(!!as.name(col_use_in_PXS))
var_tbl <- coeffs %>% left_join(fields, by="term") %>%
  filter(disease == .GlobalEnv$disease) %>%
  arrange(use_type, term)

pheno_wide <- pheno[c("FID","IID",var_tbl$term)] %>%
  drop_na() # removes NA values
#View(tibble(a = (var_tbl$term), b = (colnames(pheno_wide)[-(1:2)])) %>% mutate(c=(a==b)))
# exposures <- levels(as.factor(exposures_tbl$field.x))
# 
# # creates a second pheno table for just exposures and removes any individuals
# # with NA values for any of the exposures
# pheno_wide <- pheno %>% select(FID,IID,all_of(exposures)) %>% drop_na()
# removed_IIDs <- pheno %>% filter(!(IID %in% pheno_wide$IID)) %>% select(FID,IID)
# loc_out <- paste0(dir_out,"IIDs_NA_exposures.txt")
# write.table(removed_IIDs, loc_out, sep=" ", quote=FALSE, row.names=FALSE)
# pheno <- pheno %>% filter(IID %in% pheno_wide$IID)
# # expands the pheno table to have a column for each possible categorical answer
# # also sets negative quantitative values to 0 so they have no weight in the PXS
# for (exposure in exposures) {
#   slice <- exposures_tbl %>% filter(field.x == exposure)
#   
#   if (slice$data_type[1] == "quantitative") {
#     col_name <- slice$field.x[1]
#     column <- pheno[col_name]
#     column[column < 0 | is.na(column)] <- 0
#     pheno_wide[col_name] <- column
#     print(paste("Added column for",col_name))
#   } else if (slice$data_type[1] == "categorical") {
#     the_field <- slice$field.x[1]
#     field_codes <- slice$coding
#     for (field_code in field_codes) {
#       col_name <- paste0(the_field,"_",field_code)
#       column <- ifelse((pheno[the_field] == field_code) &
#                          (!is.na(pheno[the_field])),1,0)
#       pheno_wide[col_name] <- column
#       print(paste("Added column for",col_name))
#     }
#   }
# }

# created vector of just the XWAS exposure weights (coefficients)
#coeffs_vec <- exposures_tbl[[col_coeff]]
coeffs_vec <- var_tbl$estimate
coeffs_vec[is.na(coeffs_vec)] <- 0
num_cols_skip <- ncol(pheno_wide) - nrow(var_tbl)
PXSs <- c()
# computes dot product of expanded phenotype coding and XWAS coefficients for
# each individual in pheno file
for (i in 1:nrow(pheno_wide)) {
  iid_pheno_vec <- pheno_wide[i,(num_cols_skip+1):ncol(pheno_wide)] %>%
    unlist(use.names=FALSE)
  
  PXS <- (coeffs_vec %*% iid_pheno_vec)[1]
  PXSs <- c(PXSs,PXS)
  
  if (i %% 5000 == 0) {
    print(paste0("Computed PXS for ",i/1000,"k individuals"))
  }
}
# saves shortened version of pheno table with just IIDs and PXS
col_PXS <- paste0("PXS_",disease)
out_PXS <- pheno_wide %>% select(FID,IID)
out_PXS[col_PXS] <- PXSs

# checks normality
ggplot(out_PXS, aes(x=PXS_T2D)) +
  geom_density( color="red") +
  stat_function(fun = dnorm, n = 101, args = list(mean = mean(out_PXS$PXS_T2D), sd = sd(out_PXS$PXS_T2D)))
# forces standardized PXS value:
out_PXS[col_PXS] <- (PXSs - mean(PXSs)) / sd(PXSs)


loc_out <- paste0(dir_out,"PXS_",disease,".txt")
write.table(out_PXS,loc_out,sep=" ", row.names=FALSE, quote=FALSE)

print("Done computing PXS")

### Computes env LM: PXS ~ sex + age + assessment_center + PCs
PXS_lm_tbl <- pheno %>% select(FID,IID,all_of(col_covs)) %>%
  right_join(out_PXS, by=c("FID","IID")) %>%
  rename("PXS" = all_of(col_PXS)) %>%
  mutate(sex = as.factor(sex),
         assessment_center = as.factor(assessment_center)) %>%
  select(-FID,-IID) %>%
  drop_na()

PXS_lm <- lm(data=PXS_lm_tbl, PXS ~ .)
loc_out <- paste0("PXS_",disease,"_envLM.rds")
saveRDS(PXS_lm,loc_out)
summary(PXS_lm)

print("Done computing envLM")

### Computes Random Forest
#set.seed(1)
#pheno.rf <- randomForest(PXS ~ . , data = PXS_lm_tbl,
#                         importance=TRUE, do.trace=TRUE)
#loc_out <- paste0("PXS_",disease,"_RF.rds")
#saveRDS(pheno.rf,loc_out)
#importance(pheno.rf)
#pheno.rf
