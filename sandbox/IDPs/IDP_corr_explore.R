library(tidyverse)
library(data.table)
source('code/paths.R')
source('../ePC/code/helper_functions.R')

UKB_dict <- as_tibble(fread('../key_data/Data_Dictionary_Showcase.tsv'))

cats.IDP <- c(110, 112, 119, 106, 107, 111, 109)

fieldIDs.IDP <- (UKB_dict %>% filter(Category %in% cats.IDP))$FieldID

# Load ####
raw_phenofiles <- c("/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv",
                    "/n/no_backup2/patel/uk_biobank/main_data_34521/ukb34521.csv",
                    "/n/no_backup2/patel/uk_biobank/main_data_9512/ukb9512.csv",
                    "/n/groups/patel/uk_biobank/project_22881_678133/ukb678133.csv")
pheno.IDPs.raw <- as_tibble(fread('/n/groups/patel/uk_biobank/project_22881_678133/ukb678133.csv'))
pheno.IDPs <- pheno.IDPs.raw %>%
  select(IID=eid, ends_with('-2.0') & where(is.numeric)) %>%
  filter(rowSums(!is.na(select(., -1))) > 0)
colnames(pheno.IDPs)[-1] <- paste0('f',str_replace(colnames(pheno.IDPs)[-1], '-2.0',''))

phenoEC <- as_tibble(fread('scratch/pheno_EC.txt'))

pheno <- pheno.IDPs %>% 
  left_join(phenoEC, by='IID')

cols.IDP <- colnames(pheno.IDPs)[-1]


corrs.tbl <- tibble(field1='',field2='',r=0,r.p=0)[0,]
field1 <- 'PXS_T2D'
for (col.IDP in cols.IDP) {
  temp.tbl <- pheno[,c(field1,col.IDP)] %>% drop_na()
  cor1 <- cor.test(temp.tbl[[field1]],temp.tbl[[col.IDP]])
  
  corrs.tbl <- corrs.tbl %>%
    add_row(field1=field1,field2=col.IDP,
            r=cor1$estimate, r.p=cor1$p.value)
}
corrs.tbl$r.p.adj <- p.adjust(corrs.tbl$r.p, method = 'bonferroni')

corrs.tbl <- corrs.tbl %>%
  drop_na() %>%
  left_join(UKB_dict %>%
              mutate(field2=paste0('f',FieldID))%>%
              select(field2,Field, Category),
            by='field2') %>%
  arrange(-abs(r))

fwrite(corrs.tbl, 'scratch/IDPs/PXST2D_IDPs.corr.tbl.txt', sep='\t')
hist(corrs.tbl$r)


corrs.tbl %>% filter(Category==106) %>% View()
corrs.tbl %>% filter(str_detect(Field, "Volume of")) %>% View()
