# Used to define repeatedly-used directories and paths for other code so that
# you don't have to set them in every R file.
# No need to run independently
# Load by using source('paths.R')

dir_repo <- '/n/groups/patel/nuno/PXS-pipeline/'
dir_script <- "~/group_nuno/PXS-pipeline/code/"
dir_scratch <- "~/group_nuno/PXS-pipeline/scratch/"
dir_results <- paste0(dir_repo,"final_results/")
dir_data_showcase <- "~/group_nuno/key_data/" # contains 'Data_Dictionary_Showcase.tsv' and 'Codings.tsv' from UKBB
loc_pheno <- paste0(dir_scratch,"pheno.txt")