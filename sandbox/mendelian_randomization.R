library(tidyverse)
library(data.table)
library(MendelianRandomization)

dir_T2D <- "/n/groups/patel/nuno/PXS-pipeline/scratch/T2D/"
LMM_PXS <- as_tibble(fread(paste0(dir_T2D,"LMM_PXS_T2D_bgen.txt"))) %>%
  mutate(uniqID1 = paste(CHR,BP, ALLELE1, ALLELE0, sep=":"),
         uniqID2 = paste(CHR,BP, ALLELE0, ALLELE1, sep=":"))
LMM_T2D_onset <- as_tibble(fread(paste0(dir_T2D,"LMM_T2D_onset_bgen.txt"))) %>%
  mutate(uniqID1 = paste(CHR,BP, ALLELE1, ALLELE0, sep=":"),
         uniqID2 = paste(CHR,BP, ALLELE0, ALLELE1, sep=":"))

PXS_IndSig <- as_tibble(fread("/n/groups/patel/nuno/PXS-pipeline/scratch/FUMA_results/FUMA_job228821/IndSigSNPs.txt"))
PXS_IndSig_SNPs <- PXS_IndSig$uniqID

LMM_joint1 <- LMM_PXS %>%
  filter(uniqID1 %in% PXS_IndSig_SNPs | uniqID2 %in% PXS_IndSig_SNPs) %>%
  select(SNP, PXS_BETA = BETA, PXS_BETA_SE = SE) %>%
  left_join(
    LMM_T2D_onset %>%
      filter(uniqID1 %in% PXS_IndSig_SNPs | uniqID2 %in% PXS_IndSig_SNPs) %>%
      select(SNP, T2D_BETA = BETA, T2D_BETA_SE = SE),
    by="SNP")

MRInputObject1 <- mr_input(bx = LMM_joint1$PXS_BETA,
                           bxse = LMM_joint1$PXS_BETA_SE,
                           by = LMM_joint1$T2D_BETA,
                           byse = LMM_joint1$T2D_BETA_SE,
                           snps = LMM_joint1$SNP,
                           exposure = "PXS_T2D",
                           outcome = "T2D_onset")
MRAllObject1_all <- mr_allmethods(MRInputObject1, method = "all")
MRAllObject1_all

# T2D_onset --> PXS_T2D ####
T2D_IndSig <- as_tibble(fread("/n/groups/patel/nuno/PXS-pipeline/scratch/FUMA_results/FUMA_job265315/IndSigSNPs.txt"))
T2D_IndSig_SNPs <- T2D_IndSig$uniqID

LMM_joint2 <- LMM_PXS %>%
  filter(uniqID1 %in% T2D_IndSig_SNPs | uniqID2 %in% T2D_IndSig_SNPs) %>%
  select(SNP, PXS_BETA = BETA, PXS_BETA_SE = SE) %>%
  left_join(
    LMM_T2D_onset %>%
      filter(uniqID1 %in% T2D_IndSig_SNPs | uniqID2 %in% T2D_IndSig_SNPs) %>%
      select(SNP, T2D_BETA = BETA, T2D_BETA_SE = SE),
    by="SNP")

MRInputObject2 <- mr_input(by = LMM_joint2$PXS_BETA,
                           byse = LMM_joint2$PXS_BETA_SE,
                           bx = LMM_joint2$T2D_BETA,
                           bxse = LMM_joint2$T2D_BETA_SE,
                           snps = LMM_joint2$SNP,
                           outcome = "PXS_T2D",
                           exposure = "T2D_onset")
MRAllObject2_all <- mr_allmethods(MRInputObject2, method = "all")
MRAllObject2_all


# PXS_T2D --> T2D_all ####
LMM_T2D_all <- as_tibble(fread(paste0(dir_T2D,"LMM_T2D_all_bgen.txt"))) %>%
  mutate(uniqID1 = paste(CHR,BP, ALLELE1, ALLELE0, sep=":"),
         uniqID2 = paste(CHR,BP, ALLELE0, ALLELE1, sep=":"))

LMM_joint3 <- LMM_PXS %>%
  filter(uniqID1 %in% PXS_IndSig_SNPs | uniqID2 %in% PXS_IndSig_SNPs) %>%
  select(SNP, PXS_BETA = BETA, PXS_BETA_SE = SE) %>%
  left_join(
    LMM_T2D_all %>%
      filter(uniqID1 %in% PXS_IndSig_SNPs | uniqID2 %in% PXS_IndSig_SNPs) %>%
      select(SNP, T2D_BETA = BETA, T2D_BETA_SE = SE),
    by="SNP")

MRInputObject3 <- mr_input(bx = LMM_joint3$PXS_BETA,
                           bxse = LMM_joint3$PXS_BETA_SE,
                           by = LMM_joint3$T2D_BETA,
                           byse = LMM_joint3$T2D_BETA_SE,
                           snps = LMM_joint3$SNP,
                           exposure = "PXS_T2D",
                           outcome = "T2D_all")
MRAllObject3_all <- mr_allmethods(MRInputObject3, method = "all")
MRAllObject3_all

# T2D_onset --> PXS_T2D ####
T2Dall_IndSig <- as_tibble(fread("/n/groups/patel/nuno/PXS-pipeline/scratch/FUMA_results/FUMA_job265338/IndSigSNPs.txt"))
T2Dall_IndSig_SNPs <- T2Dall_IndSig$uniqID

LMM_joint4 <- LMM_PXS %>%
  filter(uniqID1 %in% T2Dall_IndSig_SNPs | uniqID2 %in% T2Dall_IndSig_SNPs) %>%
  select(SNP, PXS_BETA = BETA, PXS_BETA_SE = SE) %>%
  left_join(
    LMM_T2D_all %>%
      filter(uniqID1 %in% T2Dall_IndSig_SNPs | uniqID2 %in% T2Dall_IndSig_SNPs) %>%
      select(SNP, T2D_BETA = BETA, T2D_BETA_SE = SE),
    by="SNP")

MRInputObject4 <- mr_input(by = LMM_joint4$PXS_BETA,
                           byse = LMM_joint4$PXS_BETA_SE,
                           bx = LMM_joint4$T2D_BETA,
                           bxse = LMM_joint4$T2D_BETA_SE,
                           snps = LMM_joint4$SNP,
                           outcome = "PXS_T2D",
                           exposure = "T2D_all")
MRAllObject4_all <- mr_allmethods(MRInputObject4, method = "all")
MRAllObject4_all