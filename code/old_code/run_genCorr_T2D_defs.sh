# Runs genetic correlation for T2D_onset and T2D_all, different definitions of T2D

# sets directories and paths (include last forward slash / )
dir_script="/home/nur479/group_nuno/PXS-pipeline/code/"
dir_scratch="/home/nur479/group_nuno/PXS-pipeline/scratch/"

pheno="T2D"

subfolder=$(echo ${dir_scratch}${pheno})
cd ${subfolder}

# genCorr(T2D_all, T2D_onset)

#########################################
## Creates BOLT-REML script and submits##
#########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 2-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o '${pheno}'_all_onset_genCorr.out
#SBATCH -e '${pheno}'_all_onset_genCorr.err

~/bolt \
--numThreads 20 \
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_scratch}'phenoEC_fullT2D.txt \
--phenoCol '${pheno}'_all \
--phenoCol '${pheno}'_onset \
--covarFile '${dir_scratch}'phenoEC_fullT2D.txt \
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
--statsFileBgenSnps '${subfolder}'/'${pheno}'_all__bgen.txt

' > ${subfolder}/genCorr_${pheno}_all_onset.sh
sbatch ${subfolder}/genCorr_${pheno}_all_onset.sh
echo Submitted genCorr for ${pheno}_all and ${pheno}_onset
