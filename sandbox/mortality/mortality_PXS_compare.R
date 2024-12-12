library(tidyverse)
library(data.table)
library(ggrepel)

# XWAS coeffs ####
XWAS.coeffs.T2D <- as_tibble(fread('input_data/xwas_coefficients.txt'))
XWAS.coeffs.death <- as_tibble(fread('scratch/mortality/xwas_coefficients.txt'))

XWAS.coeffs.both <- XWAS.coeffs.T2D %>%
  select(z, fdr, Field, FieldID, Value, FieldName, Coding, Meaning) %>%
  inner_join(XWAS.coeffs.death %>%
              select(z, fdr, Field, FieldID, Value, FieldName, Coding, Meaning),
            by=c('Field', 'FieldID', 'Value', 'FieldName', 'Coding', 'Meaning'),
            suffix = c('.T2D','.death')) %>%
  select(starts_with('z'),starts_with('fdr'), everything()) %>%
  mutate(concordant_sign = (z.T2D * z.death) > 0,
         sig.T2D = fdr.T2D < 0.05,
         sig.death = fdr.death < 0.05,
         concordant_sig = !xor(sig.T2D,sig.death))

path.out <- 'scratch/mortality/xwas_coeffs.comparison.txt'
fwrite(XWAS.coeffs.both, path.out, sep='\t')

# all vars ####
cor1 <- cor.test(XWAS.coeffs.both$z.T2D, XWAS.coeffs.both$z.death)

ggplot(XWAS.coeffs.both,
       aes(x=z.T2D,y=z.death)) +
  geom_abline(slope=1, color='gray', linetype='dashed') +
  geom_hline(yintercept = 0, color='gray', linetype='dashed') +
  geom_vline(xintercept = 0, color='gray', linetype='dashed') +
  geom_point(aes(color=concordant_sign)) +
  geom_smooth(method='lm', formula = 'y ~ x') +
  #geom_text_repel(aes(label=FieldName)) +
  #scale_x_log10() + scale_y_log10() +
  theme_light() +
  theme(legend.position = 'none') +
  labs(title = 'XWAS Z-scores for T2D vs mortality',
       subtitle = paste0('r = ', round(cor1$estimate,3)))


# one var per FieldID ####
XWAS.coeffs.both.family <- XWAS.coeffs.both %>%
  group_by(FieldID, FieldName) %>%
  summarize(z.abs.T2D = max(abs(z.T2D)),
            z.abs.death = max(abs(z.death)))


cor2 <- cor.test(XWAS.coeffs.both.family$z.abs.T2D, XWAS.coeffs.both.family$z.abs.death)

ggplot(XWAS.coeffs.both.family,
       aes(x=z.abs.T2D,y=z.abs.death)) +
  geom_abline(slope=1, color='gray', linetype='dashed') +
  geom_hline(yintercept = 0, color='gray', linetype='dashed') +
  geom_vline(xintercept = 0, color='gray', linetype='dashed') +
  geom_point() +
  geom_smooth(method='lm', formula = 'y ~ x') +
  #geom_text_repel(aes(label=FieldName)) +
  #scale_x_log10() + scale_y_log10() +
  theme_light() +
  theme(legend.position = 'none') +
  labs(title = 'XWAS Z-scores for T2D vs mortality',
       subtitle = paste0('r = ', round(cor2$estimate,3),' :: only most significant exposure per FieldID'))


#
PXS.coeffs.death <- as_tibble(fread('scratch/mortality/PXS_coefficients.txt'))
View(PXS.coeffs.death %>% filter(substring(term,1,1)=='f'))


# PXS value itself comparison ####
pheno.death <- as_tibble(fread('scratch/mortality/pheno_EC.death.txt'))
pheno.T2D <- as_tibble(fread('scratch/phenoEC_fullT2D.txt'))

pheno.both <- pheno.death %>%
  left_join(pheno.T2D %>% select(IID, contains('T2D'), contains('TxD')), by='IID')

pheno.both.noNA <- pheno.both %>%
  filter(!is.na(PXS_death), !is.na(PXS_T2D)) %>%
  mutate(died.text = ifelse(died,'Dead','Alive'),
         T2D.text = ifelse(T2D_after_assessment==1,'T2D','No T2D'))

## full ####
cor1 <- cor.test(pheno.both.noNA$PXS_death, pheno.both.noNA$PXS_T2D)

ggplot(pheno.both, aes(x=PXS_T2D, y=PXS_death)) +
  geom_point(alpha=0.05, shape=1) +
  geom_abline(slope=1) +
  geom_smooth(method='lm',formula='y~x') +
  theme_light() +
  labs(x = 'PXS-T2D', y='PXS-death',
       title = 'Comparing PXS values for T2D and death',
       subtitle = paste0('N = ', nrow(pheno.both.noNA),
                         ' :: r = ', round(cor1$estimate, 3)))

## split ####
ggplot(pheno.both.noNA, aes(x=PXS_T2D, y=PXS_death)) +
  geom_point(alpha=0.2, shape=1) +
  geom_abline(slope=1) +
  geom_smooth(method='lm',formula='y~x') +
  facet_grid(died.text ~ T2D.text) +
  theme_light() +
  labs(x = 'PXS-T2D', y='PXS-death',
     title = 'Comparing PXS values for T2D and death, split by outcome',
     subtitle = 'Only T2D cases after first assessment shown')

pheno.both.noNA %>%
  group_by(died.text, T2D.text) %>%
  summarize(PXS_death = mean(PXS_death),
            PXS_T2D = mean(PXS_T2D),
            N = n())

fwrite(pheno.both.noNA, 'scratch/mortality/pheno_EC.death.T2D.txt', sep='\t')
