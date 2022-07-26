# sets directories and paths
loc_phenos="phenotypes.txt"
dir_out="/home/nur479/scratch3/PXS_pipeline/" #include last forward slash /
dir_script="/home/nur479/jobs/PXS_pipeline/"

lines=$(cat $loc_phenos)
for disease in $lines
do

cd ${dir_script}

subfolder=$(echo ${dir_out}${disease})

cd ${subfolder}

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
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${subfolder}'/IIDs_NA_exposures.txt \
--remove '${dir_script}'bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${subfolder}'/PXS_'${disease}'.txt \
--phenoCol PXS_'${disease}' \
--covarFile '${dir_out}'pheno_EC.txt \
--covarCol f.31.0.0 \
--covarCol f.54.0.0 \
--qCovarCol f.34.0.0 \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--lmm \
--verboseStats \
--statsFile '${subfolder}'/LMM_'${disease}'.txt \
--bgenFile /n/no_backup2/patel/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/LMM_'${disease}'_bgen.txt

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
#SBATCH -o '${disease}'_PXS_BOLTREML.out
#SBATCH -e '${disease}'_PXS_BOLTREML.err

~/bolt \
--numThreads 20 \
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${subfolder}'/IIDs_NA_exposures.txt \
--remove '${dir_script}'bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${subfolder}'/PXS_'${disease}'.txt \
--phenoCol PXS_'${disease}' \
--covarFile '${dir_out}'pheno_EC.txt \
--covarCol f.31.0.0 \
--covarCol f.54.0.0 \
--qCovarCol f.34.0.0 \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--reml \
--remlNoRefine \
--bgenFile /n/no_backup2/patel/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/LMM_'${disease}'_bgen.txt

' > ${subfolder}/${disease}_PXS_BOLTREML.sh
sbatch ${subfolder}/${disease}_PXS_BOLTREML.sh
echo 'Submitted BOLT-REML for '${disease}

done
