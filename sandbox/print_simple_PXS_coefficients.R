library(tidyverse)
library(data.table)
coeffs_raw <- as_tibble(fread("../input_data/PXS_coefficients.txt"))
fields_tbl <- as_tibble(fread("../scratch/fields_tbl.txt"))
expos <- (fields_tbl %>% filter(use_type=="exposure"))$term
coeffs <- coeffs_raw %>% 
  filter(term %in% expos) %>% rowwise() %>%
  mutate(Field_Name = paste0(ifelse(fieldname=="",term,fieldname),
                             ifelse(Meaning=="","",paste0(": ", Meaning)))) %>%
  select(Field_Name, p.value, estimate) %>%
  filter(substring(Field_Name,1,2) != "pc")

loc_out <- "../sandbox/scratch/PXS_T2D_coeffs_tbl_simple.txt"
fwrite(coeffs, loc_out, sep="\t")  
