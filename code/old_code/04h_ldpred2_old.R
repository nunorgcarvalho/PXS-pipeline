# code for running LDpred2 on BRS GWAS summary stats
# and then compute PRS

library(tidyverse)
library(data.table)
library(bigsnpr)
source('code/00_paths.R')
source('/n/scratch/users/n/nur479/snp_cor/helper_functions.R')

dir_LDmatrix <- '/n/groups/patel/nuno/UKBB_EUR_LDmatrix/LD_matrices/'
LD_file_prefix <- 'UKBB_EUR_LDmatrix_chr' # followed by chr and file type
dir_geno_files <- '/n/scratch/users/n/nur479/UKBB_geno_files/'
geno_file_prefix <- 'UKBB_geno_chr' # followed by chr and file type
(NCORES <- nb_cores()) # make sure cores have been enabled, see ?assert_cores()
dir_tmp <- '/n/scratch/users/n/nur479/tmp_ldpred2/'
dir_ldpred2 <- paste0(dir_scratch,'ldpred2/')

h2 <- 0.175 # ~ldsc's h2 estimate for BRS

dir.create(dir_tmp, showWarnings = FALSE)
dir.create(dir_ldpred2, showWarnings = FALSE)

# used to get effective sample sizes
pheno <- as_tibble(fread(paste0(dir_scratch,'cohorts/BRS_cohort_ALL.txt')))
fam <- as_tibble(fread(paste0(dir_geno_files, geno_file_prefix, 1, '.fam'))) %>%
  left_join(pheno %>% select(IID, sample_group, `BRS-ALL-cov_bvr`), by=c('V2'='IID'))
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

# QC
QC <- (sumstats %>% mutate(MAF=0, rsid='',
                           `_NUM_ID_.ss`=0,`_NUM_ID_`=0))[0,]
for (chrom in 1:22) {
  print(paste0('Chr',chrom, ' :: getting MAF for QC'))
  loc_geno <- paste0(dir_geno_files, geno_file_prefix, chrom)
  obj.bigSNP <- read_bigSNP(loc_geno)
  #G_objs[[chrom]] <- obj.bigSNP
  
  sumstats_chr <- sumstats %>% filter(chr == chrom)
  map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0")) %>% as_tibble()
  df_beta_chr <- as_tibble(snp_match(sumstats_chr, map))
  
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

chrs <- c(3:5,7:22)
#chrs <- 1:22
max_i_corr <- 0
# loop ####
for (chrom in chrs) {
  print(paste0('Chr',chrom, ' Reading summary stats'))
  loc_geno <- paste0(dir_geno_files, geno_file_prefix, chrom)
  LD_tbl_chr <- as_tibble(fread(paste0(dir_LDmatrix, LD_file_prefix, chrom, '.tsv')))
  sumstats_chr <- sumstats %>%
    inner_join(LD_tbl_chr %>% select(chr=chromosome, pos=position,
                                     a0 = allele2, a1=allele1, rsid, chrom_i0)) %>%
    inner_join(QC %>% filter(chr==chrom,G_sd_outlier==FALSE) %>%
                 select(chr,pos,a0,a1,rsid))
  M_chr <- nrow(sumstats_chr)
  bim <- as_tibble(fread(paste0(loc_geno,'.bim')))
  i_SNP <- which(bim$V2 %in% sumstats_chr$rsid)
  
  # reads genotype data
  print(paste0('Chr', chrom, ': Reading bed file'))
  if (!file.exists(paste0(loc_geno, '.rds'))) {
    snp_readBed2(bedfile = paste0(loc_geno, '.bed'), backingfile = loc_geno,
                 ind.col = i_SNP, ncores=NCORES)
  }
  obj.bigSNP <- snp_attach(paste0(loc_geno, '.rds'))
  
  map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0")) %>% as_tibble()
  df_beta_chr <- as_tibble(snp_match(sumstats_chr, map)) %>%
    mutate(i_SNP_corr = max_i_corr + chrom_i0 + 1)
  
  # reads correlation matrix data
  print(paste0('Chr',chrom, ' Reading LD matrix'))
  corr0 <- readRDS(paste0(dir_LDmatrix, LD_file_prefix, chrom, '.rds'))
  max_i_corr <- max_i_corr + ncol(corr0)
  if (chrom == chrs[1]) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, backingfile = tempfile(tmpdir=dir_tmp), compact=TRUE)
    df_beta <- df_beta_chr
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
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

  
# ldpred2-auto ####
# WARNING: I've had a lot of convergence issues so I am opting for the grid method

coef_shrink <- 0.95  # reduce this up to 0.4 if you have some (large) mismatch with the LD ref
#set.seed(1)  # to get the same result every time
print('Running ldpred2')
# takes less than 2 min with 4 cores
multi_auto <- snp_ldpred2_auto(
  corr, df_beta, h2_init = ldsc_h2_est,
  ind.corr = df_beta$i_SNP_corr,
  vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES), ncores = NCORES,
  use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
  allow_jump_sign = FALSE, shrink_corr = coef_shrink)

(h2_ests <- sapply(multi_auto, function(x) x$h2_est))
# filtering out bad chains
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

# final PGS weights
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
PGS_weights_tbl <- df_beta %>%
  #select(rsid,chr,pos,a0,a1,beta,beta_se) %>%
  mutate( PGS_weight = beta_auto )
# saves to system
loc_out <- paste0(dir_ldpred2, 'BRS_PGS_weights.tsv')
fwrite(PGS_weights_tbl, loc_out, sep='\t')

PGS_weights_tbl %>% arrange(-abs(PGS_weight)) %>%
  select(PGS_weight, everything()) %>% print(n=30)
# filter out very extreme weights
PGS_weights_tbl <- PGS_weights_tbl %>%
  filter(abs(PGS_weight) < 0.5)

#

PGS_mat <- matrix(0, nrow = length(i_test), ncol = max(chrs))
for (chrom in chrs) {
  print(chrom)
  loc_geno <- paste0(dir_geno_files, geno_file_prefix, chrom)
  obj.bigSNP <- snp_attach(paste0(loc_geno, '.rds'))
  G <- obj.bigSNP$genotypes
  G_imp <- snp_fastImputeSimple(G, method = 'random', ncores=NCORES)
  
  pred_auto <- big_prodVec(G_imp, PGS_weights_tbl$PGS_weight[PGS_weights_tbl$chr==chrom],
                           ind.row = i_test,
                           ind.col = PGS_weights_tbl[PGS_weights_tbl$chr==chrom,'_NUM_ID_'][[1]])
  PGS_mat[,chrom] <- pred_auto
  #pcor(pred_auto, fam$`BRS-ALL-cov_bvr`[i_test], NULL)
}

PGS_values <- rowSums(PGS_mat)
saveRDS(PGS_values, file=paste0(dir_ldpred2,'PGS_values.rds'))
pcor(PGS_values, fam$`BRS-ALL-cov_bvr`[i_test], NULL)
