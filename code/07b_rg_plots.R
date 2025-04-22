# makes plots related to ldsc genetic correlations (and phenotypic)

# loads directories and packages ####
library(tidyverse)
library(data.table)
source('code/00_paths.R')
source('code/00_plotting.R')

shortnames <- as_tibble(fread(paste0(dir_repo,'input_data/term_shortnames.tsv'))) %>%
  select(term, shortname, category)
cat_rle <- rle(shortnames$category)
cat_tbl <- tibble(category=cat_rle$values,
                  i_stop = cumsum(cat_rle$lengths),
                  i_start = c(1,head(i_stop + 1,-1)) )
# adds dashed lines to denote categories
add_cat_lines <- function(gg, i_start, dash_size=0.5) {
  for (i in 1:length(i_start)) {
    gg <- gg +
      annotate('segment', x=i_start[i]-0.5, xend=i_start[i]-0.5,
               y=-Inf,yend=Inf, linetype='dashed', color = "gray40", dash_size) +
      annotate('segment', y=i_start[i]-0.5, yend=i_start[i]-0.5,
               x=-Inf,xend=Inf, linetype='dashed', color = "gray40", dash_size)
  }
  return(gg)
}


# phenotypic correlation ####
col_BRS <- 'BRS-ALL-cov_bvr'
shortnames_corr <- shortnames %>% filter(term == col_BRS | category !='') %>%
  arrange(term != col_BRS) %>% select(!category)
boldings <- shortnames_corr$term == col_BRS
corr_tbl <- as_tibble(fread(paste0(dir_scratch,'general_results/behavior_corr_table.tsv'))) %>%
  left_join(shortnames %>% select(term1=term,shortname1=shortname), by='term1') %>%
  left_join(shortnames %>% select(term2=term,shortname2=shortname), by='term2') %>%
  mutate(shortname1 = factor(shortname1, levels=shortnames_corr$shortname),
         shortname2 = factor(shortname2, levels=shortnames_corr$shortname))

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
    panel.grid = element_blank(),
    axis.text = element_text(face = ifelse(boldings, 'bold','plain'))
  )
gg1 <- add_cat_lines(gg1, cat_tbl$i_start+1)
gg1



# genotypic correlations ####
rg_tbl <- as_tibble(fread(paste0(dir_scratch,'general_results/behavior_rg_table.tsv'))) %>%
  select(-starts_with('h2g'),-starts_with('intercept'), -file) %>%
  mutate(across(starts_with('term'), function(x) str_replace_all(x,'_','\\.')))
rg_tbl <- rg_tbl %>% add_row(rg_tbl %>% rename(term1=term2,term2=term1)) %>%
  mutate(across(c(term1,term2), function(x) str_replace(x, 'BRS', col_BRS)))

rg_tbl1 <- rg_tbl %>% filter(term1 %in% shortnames_corr$term, term2 %in% shortnames_corr$term) %>%
  left_join(shortnames_corr %>% select(term1=term, shortname1=shortname), by='term1') %>%
  left_join(shortnames_corr %>% select(term2=term, shortname2=shortname), by='term2') %>%
  mutate(shortname1 = factor(shortname1, levels=shortnames_corr$shortname),
         shortname2 = factor(shortname2, levels=shortnames_corr$shortname))


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
    panel.grid = element_blank(),
    axis.text = element_text(face = ifelse(boldings, 'bold','plain'))
  )
gg2 <- add_cat_lines(gg2, cat_tbl$i_start+1)
gg2


# joint corr + rg plot ####
joint_tbl <- tibble(id1=0,id2=0,shortname1='',shortname2='', type='', corr=0)[0,]
for (i in 1:length(shortnames_corr$shortname)) {
  name1 <- shortnames_corr$shortname[i]
  for (j in i:length(shortnames_corr$shortname)) {
    if (j == i) {next}
    name2 <- shortnames_corr$shortname[j]
    pheno_corr <- corr_tbl$r[corr_tbl$shortname1==name1 & corr_tbl$shortname2==name2]
    gen_corr <- rg_tbl1$rg[rg_tbl1$shortname1==name1 & rg_tbl1$shortname2==name2]
    joint_tbl <- joint_tbl %>% add_row(id1=i, id2=j,
      shortname1=c(name1,name2), shortname2=c(name2,name1),
      type=c('pheno','geno'), corr = c(pheno_corr, gen_corr)
    )

  }
}
joint_tbl <- joint_tbl %>% mutate(
  across(starts_with('shortname'),
         function(x) factor(x, levels = shortnames_corr$shortname)))

gg3 <- ggplot(joint_tbl, aes(x=shortname1, y=shortname2, fill=corr)) +
  geom_tile() +
  geom_text(aes(label = sub("(^|-)0", "\\1", sprintf("%.2f", corr))), size = 1.25) +
  scale_fill_gradient2(low="#DE1B1B",mid="white", high="#93CF1A", na.value='white',
                       midpoint = 0, limits = c(-1, 1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width=20), expand = c(0, -10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width=20), expand = c(0, 0)) +
  labs(x='Behavior', y='Behavior', fill='r') +
  coord_fixed(expand = FALSE) +
  tileplot_theme +
  theme(
    axis.text = element_text(face = ifelse(boldings, 'bold','plain')),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.grid = element_blank()
  )
gg3 <- add_cat_lines(gg3, cat_tbl$i_start+1, dash_size=0.25)
gg3
loc_fig <- paste0(dir_figs,"corr_megatable")
ggsave(paste0(loc_fig,".png"), width=180, height=180, units="mm", dpi=900)

## comparing pheno vs genetic correlations ####
joint_tbl2 <- joint_tbl %>% filter(type=='pheno') %>%
  select(id1,id2, shortname1, shortname2, pheno_corr=corr) %>%
  left_join(joint_tbl %>% filter(type=='geno') %>%
              select(id1,id2, geno_corr=corr),
            by=c('id1','id2')) %>%
  mutate(corr_diff = geno_corr - pheno_corr) %>%
  arrange(-abs(corr_diff))

fwrite(joint_tbl2, paste0(dir_scratch,'general_results/pairwise_corrs.tsv'), sep='\t')

cor_corrs <- cor.test(joint_tbl2$pheno_corr, joint_tbl2$geno_corr)
ggplot(joint_tbl2, aes(x=pheno_corr, geno_corr)) +
  geom_point(shape=1) +
  geom_abline(slope=1,linetype='dashed', color=dash_color) +
  geom_hline(yintercept=0,linetype='dashed', color=dash_color) +
  geom_vline(xintercept=0,linetype='dashed', color=dash_color) +
  #geom_smooth(method='lm',formula='y~x') +
  annotate('text', label = annotate_cor.test(cor_corrs), # throws out harmless warning message
           x=Inf, y=-Inf, vjust=-0.5, hjust=1.05, parse=FALSE, size=3) +
  scale_x_continuous(limits=c(-1,1), expand=c(0,0)) +
  scale_y_continuous(limits=c(-1,1), expand=c(0,0)) +
  labs(x='Phenotypic correlation',
       y='Genetic correlation') +
  coord_fixed() +
  theme(plot.margin = margin(r = 2.5, unit = "mm"))

loc_fig <- paste0(dir_figs,"corr_comparison")
ggsave(paste0(loc_fig,".png"), width=120, height=120, units="mm", dpi=300)
