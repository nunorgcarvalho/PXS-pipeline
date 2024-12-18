# Runs genetic correlation for both PXS_T2D and T2D itself and each of its CRFs

# sets directories and paths (include last forward slash / )
dir_script="/home/nur479/group_nuno/PXS-pipeline/code/"
dir_scratch="/home/nur479/group_nuno/PXS-pipeline/scratch/"

loc_phenolist=$(echo ${dir_script}../input_data/phenotypes.txt)
phenos=$(cat $loc_phenolist)
for pheno in $phenos
do

subfolder=$(echo ${dir_scratch}${pheno})
cd ${subfolder}

loc_CRFs=$(echo ${subfolder}/${pheno}_CRFs.txt)

CRFs=$(cat $loc_CRFs)
for CRF in $CRFs
do

# genCorr(PXS_T2D, CRFs)

#########################################
## Creates BOLT-REML script and submits##
#########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 2-11:59
#SBATCH -p medium
#SBATCH --mem=100G
#SBATCH -o PXS_'${pheno}'_'${CRF}'_genCorr.out
#SBATCH -e PXS_'${pheno}'_'${CRF}'_genCorr.err

~/bolt \
--numThreads 20 \
--bed /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_scratch}'pheno_EC.txt \
--phenoCol PXS_'${pheno}' \
--phenoCol '${CRF}' \
--covarFile '${dir_scratch}'pheno_EC.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--reml \
--remlNoRefine \
--bgenFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/PXS_'${pheno}'_'${CRF}'_bgen.txt

' > ${subfolder}/genCorr_PXS_${pheno}_${CRF}.sh
sbatch ${subfolder}/genCorr_PXS_${pheno}_${CRF}.sh
echo Submitted genCorr for PXS_${pheno} and ${CRF}


# genCorr(T2D, CRFs)

#########################################
## Creates BOLT-REML script and submits##
#########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 2-11:59
#SBATCH -p medium
#SBATCH --mem=100G
#SBATCH -o '${pheno}'_'${CRF}'_genCorr.out
#SBATCH -e '${pheno}'_'${CRF}'_genCorr.err

~/bolt \
--numThreads 20 \
--bed /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_scratch}'pheno_EC.txt \
--phenoCol '${pheno}'_onset \
--phenoCol '${CRF}' \
--covarFile '${dir_scratch}'pheno_EC.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--reml \
--remlNoRefine \
--bgenFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/'${pheno}'_'${CRF}'_bgen.txt

' > ${subfolder}/genCorr_${pheno}_${CRF}.sh
#sbatch ${subfolder}/genCorr_${pheno}_${CRF}.sh
#echo Submitted genCorr for ${pheno} and ${CRF}


# genCorr(T2D_all, CRFs)

#########################################
## Creates BOLT-REML script and submits##
#########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 2-11:59
#SBATCH -p medium
#SBATCH --mem=100G
#SBATCH -o '${pheno}'_all_'${CRF}'_genCorr.out
#SBATCH -e '${pheno}'_all_'${CRF}'_genCorr.err

~/bolt \
--numThreads 20 \
--bed /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_scratch}'phenoEC_fullT2D.txt \
--phenoCol '${pheno}'_all \
--phenoCol '${CRF}' \
--covarFile '${dir_scratch}'phenoEC_fullT2D.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--reml \
--remlNoRefine \
--bgenFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/'${pheno}'_all_'${CRF}'_bgen.txt

' > ${subfolder}/genCorr_${pheno}_all_${CRF}.sh
#sbatch ${subfolder}/genCorr_${pheno}_all_${CRF}.sh
#echo Submitted genCorr for ${pheno}_all and ${CRF}

done
done
