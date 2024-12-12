# sets directories and paths (include last forward slash / )
dir_script="/home/nur479/group_nuno/PXS-pipeline/code/"
dir_scratch="/home/nur479/group_nuno/PXS-pipeline/scratch/"

loc_phenolist=$(echo ${dir_script}../input_data/phenotypes.txt)
phenos=$(cat $loc_phenolist)
for disease in $phenos
do

cd ${dir_script}

subfolder=$(echo ${dir_scratch}${disease})

cd ${subfolder}

### PXS_T2D ###
########################################
## Creates BOLT-LMM script and submits##
########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o '${disease}'_PXS_BOLTLMM.out
#SBATCH -e '${disease}'_PXS_BOLTLMM.err

~/bolt \
--numThreads 20 \
--bed /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${subfolder}'/IIDs_NA_exposures.txt \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${subfolder}'/PXS_'${disease}'.txt \
--phenoCol PXS_'${disease}' \
--covarFile '${dir_scratch}'pheno_EC.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--lmm \
--verboseStats \
--statsFile '${subfolder}'/LMM_PXS_'${disease}'.txt \
--bgenFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/LMM_PXS_'${disease}'_bgen.txt

' > ${subfolder}/${disease}_PXS_BOLTLMM.sh
sbatch ${subfolder}/${disease}_PXS_BOLTLMM.sh
echo 'Submitted BOLT-LMM for PXS_'${disease}

### PXS_T2D BMI-adj1 ###
########################################
## Creates BOLT-LMM script and submits##
########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o '${disease}'_PXS_BMIadj_BOLTLMM.out
#SBATCH -e '${disease}'_PXS_BMIadj_BOLTLMM.err

~/bolt \
--numThreads 20 \
--bed /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${subfolder}'/IIDs_NA_exposures.txt \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${subfolder}'/PXS_'${disease}'.txt \
--phenoCol PXS_'${disease}' \
--covarFile '${dir_scratch}'pheno_EC.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--qCovarCol f21001 \
--covarMaxLevels 25 \
--lmm \
--verboseStats \
--statsFile '${subfolder}'/LMM_PXS_'${disease}'_BMIadj.txt \
--bgenFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/LMM_PXS_'${disease}'_BMIadj_bgen.txt

' > ${subfolder}/${disease}_PXS_BMIadj_BOLTLMM.sh
sbatch ${subfolder}/${disease}_PXS_BMIadj_BOLTLMM.sh
echo 'Submitted BOLT-LMM for PXS_'${disease} 'BMI-adj'

### PXS_T2D BMI-adj2 ###
########################################
## Creates BOLT-LMM script and submits##
########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o '${disease}'_PXS_BMIadj2_BOLTLMM.out
#SBATCH -e '${disease}'_PXS_BMIadj2_BOLTLMM.err

~/bolt \
--numThreads 20 \
--bed /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${subfolder}'/IIDs_NA_exposures.txt \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_scratch}'pheno_EC.txt \
--phenoCol PXS_'${disease}'_BMIadj \
--covarFile '${dir_scratch}'pheno_EC.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--lmm \
--verboseStats \
--statsFile '${subfolder}'/LMM_PXS_'${disease}'_BMIadj2.txt \
--bgenFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/LMM_PXS_'${disease}'_BMIadj2_bgen.txt

' > ${subfolder}/${disease}_PXS_BMIadj2_BOLTLMM.sh
#sbatch ${subfolder}/${disease}_PXS_BMIadj2_BOLTLMM.sh
#echo 'Submitted BOLT-LMM for PXS_'${disease} 'BMI-adj2'

### T2D_all ###
########################################
## Creates BOLT-LMM script and submits##
########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o '${disease}'_all_BOLTLMM.out
#SBATCH -e '${disease}'_all_BOLTLMM.err

~/bolt \
--numThreads 20 \
--bed /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${subfolder}'/IIDs_NA_exposures.txt \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_scratch}'phenoEC_fullT2D.txt \
--phenoCol '${disease}'_all \
--covarFile '${dir_scratch}'phenoEC_fullT2D.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--lmm \
--verboseStats \
--statsFile '${subfolder}'/LMM_'${disease}'_all.txt \
--bgenFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/LMM_'${disease}'_all_bgen.txt

' > ${subfolder}/${disease}_all_BOLTLMM.sh
sbatch ${subfolder}/${disease}_all_BOLTLMM.sh
echo 'Submitted BOLT-LMM for '${disease}'_all'

### T2D_onset ###
########################################
## Creates BOLT-LMM script and submits##
########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o '${disease}'_onset_BOLTLMM.out
#SBATCH -e '${disease}'_onset_BOLTLMM.err

~/bolt \
--numThreads 20 \
--bed /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${subfolder}'/IIDs_NA_exposures.txt \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_scratch}'pheno_EC.txt \
--phenoCol '${disease}'_onset \
--covarFile '${dir_scratch}'pheno_EC.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--lmm \
--verboseStats \
--statsFile '${subfolder}'/LMM_'${disease}'_onset.txt \
--bgenFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/LMM_'${disease}'_onset_bgen.txt

' > ${subfolder}/${disease}_onset_BOLTLMM.sh
sbatch ${subfolder}/${disease}_onset_BOLTLMM.sh
echo 'Submitted BOLT-LMM for '${disease}'_onset'

done
