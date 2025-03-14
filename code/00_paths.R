# Used to define repeatedly-used directories and paths for other code so that
# you don't have to set them in every R file.
# No need to run independently
# Load by using source('paths.R')

dir_repo <- '/n/groups/patel/nuno/PXS-pipeline/'
dir_script <- "~/group_nuno/PXS-pipeline/code/"
dir_scratch <- "~/group_nuno/PXS-pipeline/scratch/"
dir_results <- paste0(dir_repo,"final_results/")
# Warning: a lot of scripts may just refer to 'scratch/', 'code/', etc., directly

dir_data_showcase <- "~/group_nuno/key_data/" # contains 'Data_Dictionary_Showcase.tsv' and 'Codings.tsv' from UKBB
loc_40PCs <- "/n/groups/patel/nuno/key_data/UKB_40PCs_500k.txt"
# the following are paths to UKB phenotype files, see `01_generate_phenofile.R`
loc_pheno_full1 <- "/n/no_backup2/patel/uk_biobank/main_data_34521/ukb34521.tab"
loc_pheno_full2 <- "/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv"
