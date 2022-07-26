# sets directories and paths
dir_out="/home/nur479/scratch3/PXS_pipeline/" #include last forward slash /
dir_script="/home/nur479/jobs/PXS_pipeline/"

loc_phenolist=$(echo ${dir_script}phenotypes_ALL.txt)

phenos=$(cat $loc_phenolist)
for pheno in $phenos
do

subfolder=$(echo ${dir_out}${pheno})
cd ${subfolder}

loc_CRFs=$(echo ${subfolder}/${pheno}_CRFs.txt)

CRFs=$(cat $loc_CRFs)
for CRF in $CRFs
do

#########################################
## Creates BOLT-REML script and submits##
#########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o PXS_'${pheno}'_'${CRF}'_genCorr.out
#SBATCH -e PXS_'${pheno}'_'${CRF}'_genCorr.err

~/bolt \
--numThreads 20 \
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_script}'bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_out}'pheno_EC.txt \
--phenoCol PXS_'${pheno}' \
--phenoCol '${CRF}' \
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
--statsFileBgenSnps '${subfolder}'/'${pheno}'_'${CRF}'_bgen.txt

' > ${subfolder}/genCorr_PXS_${pheno}_${CRF}.sh
sbatch ${subfolder}/genCorr_PXS_${pheno}_${CRF}.sh
echo Submitted genCorr for PXS_${pheno} and ${CRF}

done
done
