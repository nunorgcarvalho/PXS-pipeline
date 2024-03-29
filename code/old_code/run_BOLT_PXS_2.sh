# sets directories and paths (include last forward slash / )
dir_script="/home/nur479/jobs/PXS_pipeline/code/"
dir_scratch="/home/nur479/scratch3/PXS_pipeline/"

loc_phenolist=$(echo ${dir_script}../input_data/phenotypes.txt)
phenos=$(cat $loc_phenolist)
for disease in $phenos
do

cd ${dir_script}

subfolder=$(echo ${dir_scratch}${disease})

cd ${subfolder}

########################################
## Creates BOLT-LMM script and submits##
########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o '${disease}'_PXS_BOLTLMM2.out
#SBATCH -e '${disease}'_PXS_BOLTLMM2.err

~/bolt \
--numThreads 20 \
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${subfolder}'/IIDs_NA_exposures.txt \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${subfolder}'/PXS_'${disease}'_2.txt \
--phenoCol PXS_'${disease}'_2 \
--covarFile '${dir_scratch}'pheno_EC.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--lmm \
--verboseStats \
--statsFile '${subfolder}'/LMM_'${disease}'_2.txt \
--bgenFile /n/no_backup2/patel/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/LMM_'${disease}'_2_bgen.txt

' > ${subfolder}/${disease}_PXS_BOLTLMM.sh
sbatch ${subfolder}/${disease}_PXS_BOLTLMM.sh
echo 'Submitted BOLT-LMM for '${disease}

#########################################
## Creates BOLT-REML script and submits##
#########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o '${disease}'_PXS_BOLTREML2.out
#SBATCH -e '${disease}'_PXS_BOLTREML2.err

~/bolt \
--numThreads 20 \
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${subfolder}'/IIDs_NA_exposures.txt \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${subfolder}'/PXS_'${disease}'_2.txt \
--phenoCol PXS_'${disease}'_2 \
--covarFile '${dir_scratch}'pheno_EC.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--reml \
--remlNoRefine \
--bgenFile /n/no_backup2/patel/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/LMM_'${disease}'_2_bgen.txt

' > ${subfolder}/${disease}_PXS_BOLTREML.sh
sbatch ${subfolder}/${disease}_PXS_BOLTREML.sh
echo 'Submitted BOLT-REML for '${disease}

done
