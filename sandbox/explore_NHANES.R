library(tidyverse)
library(DBI)
library(dbplyr)
library(RSQLite)

path_to_dbname <- "scratch/NHANES/pe_summary_stats_02_2024-v2.sqlite" ## change this to your path
## connect:
con <- DBI::dbConnect(RSQLite::SQLite(), dbname=path_to_dbname)
## list tables:
DBI::dbListTables(con)

## all the variable names:
var_names <- tbl(con, "variable_names_epcf")
## all the table names (aka Data.File.Name)
table_names <- tbl(con, "table_names_epcf")
table_names


## get sample table(s) for spirometry from the NHANES 2009-2010 (F)
demo_f <- tbl(con, "DEMO_F")
spx_f <- tbl(con, "SPX_F")

spx_demo <- demo_f |> left_join(spx_f, by="SEQN")
spx_demo

spx_demo |> group_by(RIAGENDR) |> summarize(mean_fev1 = mean(SPXNFEV1))

p <- ggplot(spx_demo, aes(RIDAGEYR, SPXNFEV1, color=factor(RIAGENDR)))
p <- p + geom_point() + xlab("Age") + ylab("N_FEV1")
p




############

library(foreign)
download.file("https://wwwn.cdc.gov/nchs/nhanes/2017-2018/P_DEMO.XPT", tf <- tempfile(), mode="wb")
DEMO_I3 <- as_tibble(foreign::read.xport(tf))


############
load("/n/groups/patel/nuno/PXS-pipeline/scratch/NHANES/nhanes2018_vy.RData")

cols.merge_ID <- c('SEQN','SEQN_new','SDDSRVYR')
names.tables <- c('demographics','questionnaire','occupation','dietary',
                  'medications','chemicals','response','mortality','weights','comments')
names.tables_clean <- paste0(names.tables,'_clean')

# list.tables <- list()
# for (name.table in names.tables_clean) {
#   list.tables[[name.table]] <- as_tibble(get(name.table))
# }
data <- as_tibble(get(names.tables_clean[1]))
for (name.table in names.tables_clean[-1]){
  print(name.table)
  data <- data %>% left_join(get(name.table),
                             by = cols.merge_ID)
}
