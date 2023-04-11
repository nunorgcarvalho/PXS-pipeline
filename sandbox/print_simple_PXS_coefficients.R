coeffs_raw <- as_tibble(fread("~/jobs/PXS_pipeline/input_data/PXS_coefficients.txt"))
coeffs <- coeffs_raw %>% rowwise() %>%
  mutate(Field_Name = paste0(ifelse(is.na(fieldname),term,fieldname),
                             ifelse(is.na(Meaning),"",paste0(": ", Meaning)))) %>%
  select(Field_Name, p.value, estimate) %>%
  filter(substring(Field_Name,1,2) != "pc")

  