
library(tidyverse)
library(data.table)
library(bigsnpr)
source('code/00_paths.R')

dir_LDmatrix <- '/n/scratch/users/n/nur479/UKBB_EUR_corrmatrix/'
dir_geno_files <- '/n/scratch/users/n/nur479/UKBB_geno_files/'
geno_file_prefix <- 'UKBB_geno_chr'
pheno <- as_tibble(fread(paste0(dir_scratch,'cohorts/BRS_cohort_ALL.txt')))
N <- 10000
set.seed(2025)
IID_keep <- sample(pheno$IID, N, replace=FALSE)
NCORES <- nb_cores()
# make sure cores have been enabled, see ?assert_cores()



for (chr in 1:22) {
  loc_geno <- paste0(dir_geno_files,geno_file_prefix, chr)
  
  print(paste0('Chr', chr,': Reading bed file'))
  if (!file.exists(paste0(loc_geno,'.rds'))) {
    snp_readBed2(bedfile = paste0(loc_geno,'.bed'), backingfile = loc_geno)
  }
  obj.bigSNP <- snp_attach(paste0(loc_geno,'.rds'))
  i_keep <- which(obj.bigSNP$fam$sample.ID %in% IID_keep)
  
  G   <- obj.bigSNP$genotypes
  POS <- obj.bigSNP$map$physical.pos
  CHR <- obj.bigSNP$map$chromosome
  # computes genetic distance in centiMorgans (cM)
  POS2 <- snp_asGeneticPos(CHR, POS, dir = dir_geno_files, ncores = NCORES)
  
  print(paste0('Chr', chr,': Computing correlation matrix'))
  t1 <- Sys.time()
  # computes within 1cM window, as recommended by LDsc:
  # https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
  corr0 <- snp_cor(G, size = 1, infos.pos = POS2, ncores=NCORES,
                   ind.row = i_keep)
  print(Sys.time() - t1)
  
  print(paste0('Chr', chr,': Saving correlation matrix'))
  loc_out <- paste0(dir_LDmatrix,'UKBB_EUR_corrmatrix_chr',chr)
  saveRDS(corr0, file=paste0(loc_out,'.rds'))
  
}