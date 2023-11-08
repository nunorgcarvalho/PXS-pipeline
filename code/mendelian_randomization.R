# Libraries and paths ####
library(tidyverse)
library(data.table)
library(MendelianRandomization)

setwd("/n/groups/patel/nuno/PXS-pipeline/")
dir_T2D <- "/n/groups/patel/nuno/PXS-pipeline/scratch/T2D/"

# Shared Code ####
LMM_PXS <- as_tibble(fread(paste0(dir_T2D,"LMM_PXS_T2D_bgen.txt"))) %>%
  mutate(uniqID1 = paste(CHR,BP, ALLELE1, ALLELE0, sep=":"),
         uniqID2 = paste(CHR,BP, ALLELE0, ALLELE1, sep=":"))
LMM_T2D_onset <- as_tibble(fread(paste0(dir_T2D,"LMM_T2D_onset_bgen.txt"))) %>%
  mutate(uniqID1 = paste(CHR,BP, ALLELE1, ALLELE0, sep=":"),
         uniqID2 = paste(CHR,BP, ALLELE0, ALLELE1, sep=":"))

PXS_IndSig <- as_tibble(fread("/n/groups/patel/nuno/PXS-pipeline/scratch/FUMA_results/FUMA_job228821/IndSigSNPs.txt"))
PXS_IndSig_SNPs <- PXS_IndSig$uniqID

# PXS_T2D --> T2D onset ####

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

# PXS_T2D --> Finngen_T2D ####
library(TwoSampleMR)
avail <- available_outcomes()

Finngen_pheno <- "finn-b-T2D_WIDE"
Finngen_sf <- as_tibble(extract_outcome_data(PXS_IndSig$rsID, Finngen_pheno))
Finngen_sig <- as_tibble(extract_instruments(Finngen_pheno))

LMM_joint5 <- LMM_PXS %>%
  filter(uniqID1 %in% PXS_IndSig_SNPs | uniqID2 %in% PXS_IndSig_SNPs) %>%
  select(SNP, PXS_BETA = BETA, PXS_BETA_SE = SE) %>%
  left_join(
    Finngen_sf %>%
      select(SNP, T2D_BETA = beta.outcome, T2D_BETA_SE = se.outcome),
    by="SNP") %>%
  drop_na()

MRInputObject5 <- mr_input(bx = LMM_joint5$PXS_BETA,
                           bxse = LMM_joint5$PXS_BETA_SE,
                           by = LMM_joint5$T2D_BETA,
                           byse = LMM_joint5$T2D_BETA_SE,
                           snps = LMM_joint5$SNP,
                           exposure = "PXS_T2D",
                           outcome = "T2D_Finngen")
MRAllObject5_all <- mr_allmethods(MRInputObject5, method = "all")
MRAllObject5_all

# Fingenn_T2D --> PXS_T2D ####
LMM_joint6 <- LMM_PXS %>%
  filter(SNP %in% Finngen_sig$SNP) %>%
  select(SNP, PXS_BETA = BETA, PXS_BETA_SE = SE) %>%
  left_join(
    Finngen_sig %>%
      select(SNP, T2D_BETA = beta.exposure, T2D_BETA_SE = se.exposure),
    by="SNP") %>%
  drop_na()

MRInputObject6 <- mr_input(by = LMM_joint6$PXS_BETA,
                           byse = LMM_joint6$PXS_BETA_SE,
                           bx = LMM_joint6$T2D_BETA,
                           bxse = LMM_joint6$T2D_BETA_SE,
                           snps = LMM_joint6$SNP,
                           outcome = "PXS_T2D",
                           exposure = "T2D_Finngen")
MRAllObject6_all <- mr_allmethods(MRInputObject6, method = "all")
MRAllObject6_all


# PXS_T2D --> GIANT_BMI ####

GIANT_pheno <- "ieu-b-40"
GIANT_sf <- as_tibble(extract_outcome_data(PXS_IndSig$rsID, GIANT_pheno))
GIANT_sig <- as_tibble(extract_instruments(GIANT_pheno))

LMM_joint7 <- LMM_PXS %>%
  filter(uniqID1 %in% PXS_IndSig_SNPs | uniqID2 %in% PXS_IndSig_SNPs) %>%
  select(SNP, PXS_BETA = BETA, PXS_BETA_SE = SE) %>%
  left_join(
    GIANT_sf %>%
      select(SNP, BMI_BETA = beta.outcome, BMI_BETA_SE = se.outcome),
    by="SNP") %>%
  drop_na()

MRInputObject7 <- mr_input(bx = LMM_joint7$PXS_BETA,
                           bxse = LMM_joint7$PXS_BETA_SE,
                           by = LMM_joint7$BMI_BETA,
                           byse = LMM_joint7$BMI_BETA_SE,
                           snps = LMM_joint7$SNP,
                           exposure = "PXS_T2D",
                           outcome = "BMI_GIANT")
MRAllObject7_all <- mr_allmethods(MRInputObject7, method = "all")
MRAllObject7_all

# GIANT_BMI --> PXS_T2D ####
LMM_joint8 <- LMM_PXS %>%
  filter(SNP %in% GIANT_sig$SNP) %>%
  select(SNP, PXS_BETA = BETA, PXS_BETA_SE = SE) %>%
  left_join(
    GIANT_sig %>%
      select(SNP, BMI_BETA = beta.exposure, BMI_BETA_SE = se.exposure),
    by="SNP") %>%
  drop_na()

MRInputObject8 <- mr_input(by = LMM_joint8$PXS_BETA,
                           byse = LMM_joint8$PXS_BETA_SE,
                           bx = LMM_joint8$BMI_BETA,
                           bxse = LMM_joint8$BMI_BETA_SE,
                           snps = LMM_joint8$SNP,
                           outcome = "PXS_T2D",
                           exposure = "BMI_GIANT")
MRAllObject8_all <- mr_allmethods(MRInputObject8, method = "all")
MRAllObject8_all

# GIANT_BMI --> Finngen_T2D ####
LMM_joint9 <- GIANT_sig %>%
  select(SNP, BMI_BETA = beta.exposure, BMI_BETA_SE = se.exposure) %>%
  left_join(as_tibble(extract_outcome_data(GIANT_sig$SNP, Finngen_pheno)) %>%
              select(SNP, T2D_BETA = beta.outcome, T2D_BETA_SE = se.outcome),
            by = "SNP") %>%
  drop_na()

MRInputObject9 <- mr_input(bx = LMM_joint9$BMI_BETA,
                           bxse = LMM_joint9$BMI_BETA_SE,
                           by = LMM_joint9$T2D_BETA,
                           byse = LMM_joint9$T2D_BETA_SE,
                           snps = LMM_joint9$SNP,
                           exposure = "BMI_GIANT",
                           outcome = "Finngenn_T2D")
MRAllObject9_all <- mr_allmethods(MRInputObject9, method = "ivw")
MRAllObject9_all
