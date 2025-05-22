# code for running LDpred2 on BRS GWAS summary stats
# and then compute PRS

library(tidyverse)
library(data.table)
library(bigsnpr)
source('code/00_paths.R')
source('/n/scratch/users/n/nur479/snp_cor/helper_functions.R')

dir_LDmatrix <- '/n/scratch/users/n/nur479/snp_cor/matrices/'
LD_file_prefix <- 'UKBB_EUR_LDmatrix_chr' # followed by chr and file type
dir_geno_files <- '/n/scratch/users/n/nur479/UKBB_geno_files/'
geno_file_prefix <- 'UKBB_geno_chr' # followed by chr and file type
(NCORES <- nb_cores()) # make sure cores have been enabled, see ?assert_cores()
dir_tmp <- '/n/scratch/users/n/nur479/tmp_ldpred2/'
dir_ldpred2 <- paste0(dir_scratch,'ldpred2/')

dir.create(dir_tmp, showWarnings = FALSE)
dir.create(dir_ldpred2, showWarnings = FALSE)

# used to get effective sample sizes
pheno <- as_tibble(fread(paste0(dir_scratch,'cohorts/BRS_cohort_ALL.txt')))
fam <- as_tibble(fread(paste0(dir_geno_files, geno_file_prefix, 1, '.fam'))) %>%
  left_join(pheno %>%
              select(IID, sample_group, BRS = `BRS-ALL-cov_bvr`,
                     sex, age, starts_with('pc')), by=c('V2'='IID'))
i_train <- which(fam$sample_group=='A')
i_test <- which(fam$sample_group=='B')

# reads summary statistics
sumstats <- as_tibble(bigreadr::fread2(
  paste0(dir_scratch,'LMM/LMM.ALL_train.BRS-ALL-cov_bvr.txt'))) %>%
  select(chr=CHR,pos=BP, a0=ALLELE0,a1=ALLELE1, beta=BETA, beta_se=SE) %>%
  filter(!is.nan(beta_se), !is.nan(beta))
M_total <- nrow(sumstats)
N <- sum(pheno$sample_group=='A')
sumstats$n_eff <- N


#G_objs <- list()

# QC ####
QC <- (sumstats %>% mutate(MAF=0, rsid='',
                           `_NUM_ID_.ss`=0,`_NUM_ID_`=0))[0,]
for (chrom in 1:22) {
  print(paste0('Chr',chrom, ' :: Reading genotypes'))
  loc_geno <- paste0(dir_geno_files, geno_file_prefix, chrom)
  obj.bigSNP <- read_bigSNP(loc_geno)
  #G_objs[[chrom]] <- obj.bigSNP
  
  sumstats_chr <- sumstats %>% filter(chr == chrom)
  map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0")) %>% as_tibble()
  df_beta_chr <- as_tibble(snp_match(sumstats_chr, map))
  print(paste0('Chr',chrom, ' :: Getting MAF for QC'))
  MAF <- snp_MAF(obj.bigSNP$genotypes, ind.row = i_train, ncores=NCORES)
  QC_chr <- df_beta_chr %>% mutate(
    MAF = MAF[df_beta_chr[['_NUM_ID_']]]
  )
  QC <- QC %>% add_row(QC_chr)
}
# from ldpred2 tutorial
sd_y <- with(QC, sqrt(0.5* (n_eff * beta_se**2 + beta**2))) %>% quantile(0.01)
QC <- QC %>% mutate(
  G_sd_est = sd_y / sqrt(n_eff * beta_se**2 + beta**2),
  G_sd_obs = sqrt(2 * MAF * (1 - MAF))
)
#with(QC, cor.test(G_sd_est, G_sd_obs))
(lm1 <- lm(G_sd_obs ~ G_sd_est, data=QC)) # should be near 1
QC <- QC %>% mutate(
  G_sd_resid = lm1$residuals,
  G_sd_resid_student = G_sd_resid / sd(G_sd_resid),
  G_sd_outlier = abs(G_sd_resid) > 0.03  # following code from ldpred2 author:
) #https://github.com/privefl/paper-misspec/blob/1f459fd0d6aac66f0f2945446bfaa435b585747b/code/prepare-sumstats-finngen/T2D.R#L73C1-L74C1
(mean(QC$G_sd_outlier))
saveRDS(QC, file=paste0(dir_ldpred2,'QC.rds'))


sumstats_QC <- QC %>% filter(G_sd_outlier == FALSE) %>%
  select(chr, pos, a0, a1, beta, beta_se, n_eff)


# combine LD matrices ####
max_i_corr <- 0
chrs <- 1:22
# loop ####
for (chrom in chrs) {
  print(paste0('Chr',chrom, ' :: Reading summary stats'))
  loc_geno <- paste0(dir_geno_files, geno_file_prefix, chrom)
  sumstats_chr <- sumstats_QC %>% filter(chr == chrom)
  
  bim <- read_bim(loc_geno) %>%
    select(chr=CHR,rsid=SNP,pos=BP,a1=A1,a0=A2)
  M_chr <- nrow(bim)
  
  df_beta_chr <- as_tibble(snp_match(sumstats_chr, bim)) %>%
    mutate(i_SNP_corr = max_i_corr + `_NUM_ID_`)
  
  # reads correlation matrix data
  print(paste0('Chr',chrom, ' Reading LD matrix'))
  corr0 <- readRDS(paste0(dir_LDmatrix, LD_file_prefix, chrom, '.rds'))
  if (nrow(corr0) != M_chr) {
    print(paste0('!WARNING! Chr',chrom,' :: Number of variants in corr matrix (',
                 nrow(corr0),') does not match number of variants in bim file (',
                 M_chr,')'))
  }
  max_i_corr <- max_i_corr + M_chr
  
  print(paste0('Chr',chrom, ' Combining LD matrices'))
  if (chrom == chrs[1]) {
    ld <- Matrix::colSums(corr0^2, na.rm = TRUE)
    corr <- as_SFBM(corr0, backingfile = tempfile(tmpdir=dir_tmp), compact=TRUE)
    df_beta <- df_beta_chr
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2,, na.rm = TRUE))
    corr <- corr$add_columns(corr0, nrow(corr))
    df_beta <- df_beta %>% add_row(df_beta_chr)
  }
}

# ldsc h2 estimate
# Estimate of h2 from LD Score regression
ld_keep <- ld[df_beta$i_SNP_corr]
(ldsc <- with(df_beta, snp_ldsc(ld_keep, length(ld_keep),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))
(ldsc_h2_est <- ldsc[["h2"]])



coef_shrink <- 0.95  # reduce this up to 0.4 if you have some (large) mismatch with the LD ref
#set.seed(1)  # to get the same result every time
print('Running ldpred2')
# takes less than 2 min with 4 cores
multi_auto <- snp_ldpred2_auto(
  corr, df_beta, h2_init = ldsc_h2_est,
  ind.corr = df_beta$i_SNP_corr,
  vec_p_init = seq_log(1e-4, 0.4, length.out = NCORES), ncores = NCORES,
  #use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
  allow_jump_sign = FALSE, shrink_corr = coef_shrink)

(h2_ests <- sapply(multi_auto, function(x) x$h2_est))
# filtering out bad chains
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

# final PGS weights
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
PGS_weights_tbl <- df_beta %>%
  mutate( PGS_weight = beta_auto ) %>%
  select(rsid,chr,pos,a0,a1,beta,beta_se, PGS_weight, everything())
# saves to system
loc_out <- paste0(dir_ldpred2, 'BRS_PGS_weights.tsv')
fwrite(PGS_weights_tbl, loc_out, sep='\t')

PGS_weights_tbl %>% arrange(-abs(PGS_weight)) %>%
  select(PGS_weight, everything()) %>% print(n=10)

# computing PGS values ####
PGS_mat <- matrix(as.numeric(NA), nrow = nrow(fam), ncol = max(chrs))
for (chrom in chrs) {
  print(paste0('Chr', chrom, ' :: Computing PGS'))
  loc_geno <- paste0(dir_geno_files, geno_file_prefix, chrom)
  obj.bigSNP <- read_bigSNP(loc_geno)
  
  pred_auto <- big_prodVec(obj.bigSNP$genotypes,
                           PGS_weights_tbl$PGS_weight[PGS_weights_tbl$chr==chrom],
                           #ind.row = i_test, # compute for just test
                           ind.col = PGS_weights_tbl[['_NUM_ID_']][PGS_weights_tbl$chr==chrom])
  PGS_mat[,chrom] <- pred_auto
}

PGS_values <- rowSums(PGS_mat)
saveRDS(PGS_values, file=paste0(dir_ldpred2,'PGS_values.rds'))

# saves PGS values to pheno
fam$BRSPGS <- PGS_values
pheno <- pheno %>% left_join(fam %>% select(IID=V2, BRSPGS), by='IID')
# saves to system
loc_out <- paste0(dir_scratch,'cohorts/BRS_cohort_ALL.txt')
fwrite(pheno, loc_out, sep='\t')

# computes partial correlation (and R2) of PGS on BRS
test_tbl <- pheno %>% filter(sample_group == 'B')
col_covs <- c('sex','age',paste0('pc',1:40))
(pcor_BRSPGS_BRS <- pcor(test_tbl$BRSPGS, test_tbl$`BRS-ALL-cov_bvr`, test_tbl[,col_covs]))
(pcor_BRSPGS_BRS[1]^2)

# R2 in training group = 35.9% (!!! overfit, obviously)
# R2 in testing group = 6.26%

# R2 for CRS
(pcor_BRSPGS_CRS <- pcor(test_tbl$BRSPGS, test_tbl$`CRS-ALL`, test_tbl[,col_covs]))
(pcor_BRSPGS_CRS[1]^2)

# R2 for PRS-T2D
(pcor_BRSPGS_PRST2D <- pcor(test_tbl$BRSPGS, test_tbl$`PRS_T2D`, test_tbl[,col_covs]))
(pcor_BRSPGS_PRST2D[1]^2)
