#!/bin/sh
#SBATCH -c 20
#SBATCH -t 3-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o rg.PXS_death.PXS_T2D.out
#SBATCH -e rg.PXS_death.PXS_T2D.err

~/bolt \
--numThreads 20 \
--bed /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove /home/nur479/group_nuno/PXS-pipeline/scratch/mortality/IIDs_NA_exposures.txt \
--remove /home/nur479/group_nuno/PXS-pipeline/code/../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile /home/nur479/group_nuno/PXS-pipeline/scratch/mortality/pheno_EC.death.T2D.txt \
--phenoCol PXS_death \
--phenoCol PXS_T2D \
--covarFile /home/nur479/group_nuno/PXS-pipeline/scratch/mortality/pheno_EC.death.T2D.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol yob \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--reml \
--remlNoRefine \
--bgenFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps /home/nur479/group_nuno/PXS-pipeline/scratch/mortality/LMM.PXS_death.bgen.txt


