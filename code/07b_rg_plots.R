# makes plots related to ldsc genetic correlations (and phenotypic)

# loads directories and packages ####
library(tidyverse)
library(data.table)
source('code/00_paths.R')

shortnames <- as_tibble(fread(paste0(dir_repo,'input_data/term_shortnames.tsv'))) %>%
  select(term, shortname, category) %>%
  filter(!(term %in% c('BRS-ALL-cov_bvr','CRS-ALL','PRS_T2D')))

cat_rle <- rle(shortnames$category)
cat_tbl <- tibble(category=cat_rle$values,
                  i_stop = cumsum(cat_rle$lengths),
                  i_start = c(1,head(i_stop + 1,-1)) )

# shared theme for multiple plots
tileplot_theme <- theme(
  legend.key.size = unit(2, "mm"),
  legend.title = element_text(size=6),
  legend.text = element_text(size=5),
  legend.margin = margin(0,-1,0,-2, unit="mm"),
  legend.text.align = 1,
  plot.title = element_text(size=7),
  plot.subtitle = element_text(size=5),
  axis.title = element_text(size=6),
  axis.text = element_text(size=5)
)

# adds dashed lines to denote categories
add_cat_lines <- function(gg, i_start) {
  for (i in 1:length(i_start)) {
    gg <- gg +
      annotate('segment', x=i_start[i]-0.5, xend=i_start[i]-0.5,
               y=-Inf,yend=Inf, linetype='dashed', color = "gray40", size=0.5) +
      annotate('segment', y=i_start[i]-0.5, yend=i_start[i]-0.5,
               x=-Inf,xend=Inf, linetype='dashed', color = "gray40", size=0.5)
  }
  return(gg)
}


# phenotypic correlation ####
corr_tbl <- as_tibble(fread(paste0(dir_scratch,'general_results/behavior_corr_table.tsv'))) %>%
  left_join(shortnames %>% rename(shortname1=shortname), by=c('term1'='term')) %>%
  left_join(shortnames %>% rename(shortname2=shortname), by=c('term2'='term')) %>%
  mutate(shortname1 = factor(shortname1, levels=shortnames$shortname),
         shortname2 = factor(shortname2, levels=shortnames$shortname))

gg1 <- ggplot(corr_tbl, aes(x=shortname1, y=shortname2, fill=r)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", r)), size=2) +
  scale_fill_gradient2(low="#DE1B1B",mid="white", high="#93CF1A",
                       midpoint = 0, limits = c(-1, 1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width=20), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width=20), expand = c(0, 0)) +
  xlab('Behavior') + ylab('Behavior') +
  theme_bw() +
  tileplot_theme +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid = element_blank()
  )
gg1 <- add_cat_lines(gg1, cat_tbl$i_start)
gg1



# genotypic correlations ####
shortnames_rg <- shortnames[c(49,1:48),] # manual code, I'm sorry
terms_keep <- shortnames_rg$term
rg_tbl <- as_tibble(fread(paste0(dir_scratch,'general_results/behavior_rg_table.tsv'))) %>%
  select(-starts_with('h2g'),-starts_with('intercept'), -file) %>%
  mutate(across(starts_with('term'), function(x) str_replace_all(x,'_','\\.')))
rg_tbl <- rg_tbl %>% add_row(rg_tbl %>% rename(term1=term2,term2=term1))

rg_tbl1 <- rg_tbl %>% filter(term1 %in% terms_keep, term2 %in% terms_keep) %>%
  left_join(shortnames_rg %>% select(term1=term, shortname1=shortname), by='term1') %>%
  left_join(shortnames_rg %>% select(term2=term, shortname2=shortname), by='term2') %>%
  mutate(shortname1 = factor(shortname1, levels=shortnames_rg$shortname),
         shortname2 = factor(shortname2, levels=shortnames_rg$shortname))


gg2 <- ggplot(rg_tbl1, aes(x=shortname1, y=shortname2, fill=rg)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", rg)), size=2) +
  scale_fill_gradient2(low="#DE1B1B",mid="white", high="#93CF1A", na.value='white',
                       midpoint = 0, limits = c(-1, 1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width=20), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width=20), expand = c(0, 0)) +
  xlab('Behavior') + ylab('Behavior') +
  theme_bw() +
  tileplot_theme +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid = element_blank()
  )
gg2 <- add_cat_lines(gg2, cat_tbl$i_start+1)
gg2
