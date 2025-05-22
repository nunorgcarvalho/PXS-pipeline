# computes LD matrices. Call as an Rscript, example:
# Rscript compute_LD_chr.R 21
# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !grepl("^[0-9]+$", args[1])) {
  stop("Usage: Rscript compute_LD_chr.R <chr>\n  where <chr> is a chromosome number (1-22)")
}
chr <- as.integer(args[1])
if (!(chr %in% 1:22)) stop("Chromosome must be between 1 and 22.")

# Load paths and set parameters
library(tidyverse)
library(data.table)
library(bigsnpr)
source('code/00_paths.R')

dir_LDmatrix <- '/n/scratch/users/n/nur479/UKBB_EUR_corrmatrix/'
dir_geno_files <- '/n/scratch/users/n/nur479/UKBB_geno_files/'
geno_file_prefix <- 'UKBB_geno_chr'

pheno <- as_tibble(fread(paste0(dir_scratch, 'cohorts/BRS_cohort_ALL.txt')))
N <- 500
set.seed(2025)
IID_keep <- sample(pheno$IID, N, replace = FALSE)

NCORES <- nb_cores()

# Process chromosome
loc_geno <- paste0(dir_geno_files, geno_file_prefix, chr)

message(paste0('Chr', chr, ': Reading bed file'))
if (!file.exists(paste0(loc_geno, '.rds'))) {
  snp_readBed2(bedfile = paste0(loc_geno, '.bed'), backingfile = loc_geno)
}
obj.bigSNP <- snp_attach(paste0(loc_geno, '.rds'))
i_keep <- which(obj.bigSNP$fam$sample.ID %in% IID_keep)

G   <- obj.bigSNP$genotypes
POS <- obj.bigSNP$map$physical.pos
CHR <- obj.bigSNP$map$chromosome
POS2 <- snp_asGeneticPos(CHR, POS, dir = dir_geno_files, ncores = NCORES)

message(paste0('Chr', chr, ': Computing correlation matrix'))
t1 <- Sys.time()
corr0 <- snp_cor(G, size = 1, infos.pos = POS2, ncores = NCORES, ind.row = i_keep)
message(Sys.time() - t1)

message(paste0('Chr', chr, ': Saving correlation matrix'))
loc_out <- paste0(dir_LDmatrix, 'UKBB_EUR_corrmatrix_chr', chr)
saveRDS(corr0, file = paste0(loc_out, '.rds'))

# THE BELOW CODE GENERATES BASH SCRIPT FILES TO CALL THIS SCRIPT
# RUN THIS MANUALLY. WON'T GET RUN WHEN CALLED

if (FALSE) {
  source('code/00_paths.R')
  
  dir_LDmatrix <- '/n/scratch/users/n/nur479/UKBB_EUR_corrmatrix/'
  dir_scripts <- paste0(dir_LDmatrix, 'scripts/')
  dir.create(dir_scripts, showWarnings = FALSE)
  dir_geno_files <- '/n/scratch/users/n/nur479/UKBB_geno_files/'
  geno_file_prefix <- 'UKBB_geno_chr'
  loc_Rscript <- paste0(dir_repo,'code/03e_LDmatrix.R')
  N <- 5000
  
  scripts <- c()
  for (chr in 1:22) {
    print(paste0('Generating script for chr ', chr))
    sbatch_oe <- paste0('chr',chr,'_LDmatrix')
    loc_bim <- paste0(dir_geno_files, geno_file_prefix, chr, '.bim')
    loc_bed <- paste0(dir_geno_files, geno_file_prefix, chr, '.bed')
    M <- as.integer(system(paste0("wc -l < ",loc_bim), intern = TRUE))
    # Assumes loc_bed is a path to the file
    file_size_gb <- as.numeric(system(paste("stat -c %s", loc_bed), intern = TRUE)) / (1024^3)
    GB_RAM <- ceiling(file_size_gb) * 4 # just to be extra safe
    
    script_file <- paste0('chr',chr,'_script.sh')
    loc_script <- paste0(dir_scripts, script_file)
    scripts <- c(scripts, script_file)
    
    # for reference, w/ 20 cores, N=10k, M=11k took 12mins
    # time scales linearly w/ N, quadratically w/ M
    hours <- ceiling(6.5*(M / 11000)^2 * (N / 10000) * (12/60)) + 1
    partition <- 'short'
    if (hours > 12) {partition <- 'medium'}
    
    script <- cat('#!/bin/sh
#SBATCH -c 20
#SBATCH -t 0-',hours,':00
#SBATCH -p ',partition,'
#SBATCH --mem=',GB_RAM,'G
#SBATCH -o ',sbatch_oe,'.out
#SBATCH -e ',sbatch_oe,'.err

module load gcc/9.2.0 R/4.1.2
cd ',dir_repo,'
Rscript ',loc_Rscript,' ',chr,'
',sep='', file=loc_script)
  }
  # run from inside dir_scripts
  cat('for file in',paste0(scripts, collapse=' '),'; do sbatch "$file";done')
}

