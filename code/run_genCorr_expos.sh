# sets directories and paths (include last forward slash / )
dir_script="/home/nur479/jobs/PXS_pipeline/code/"
dir_scratch="/home/nur479/scratch3/PXS_pipeline/"

loc_phenolist=$(echo ${dir_script}../input_data/phenotypes.txt)
phenos=$(cat $loc_phenolist)
for pheno in $phenos
do

subfolder=$(echo ${dir_scratch}${pheno})
cd ${subfolder}

loc_CRFs=$(echo ${subfolder}/${pheno}_CRFs.txt)
loc_expos=$(echo ${subfolder}/${pheno}_exposures.txt)

CRFs=$(cat $loc_CRFs)
expos=$(cat $loc_expos)

for CRF in $CRFs
do

for expo in $expos
do

expo_folder=$(echo ${dir_scratch}exposures/${expo})
mkdir -p ${expo_folder}
cd ${expo_folder}

# skips if already calculated
out_file=${expo_folder}/${CRF}'_'${expo}_genCorr.out
if [[ -f "$out_file" ]]; then
echo Skipping ${CRF} and ${expo} since results already exist in folder 
continue
fi

#########################################
## Creates BOLT-REML script and submits##
#########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 2-11:59
#SBATCH -p medium
#SBATCH --mem=100G
#SBATCH -o '${CRF}'_'${expo}'_genCorr.out
#SBATCH -e '${CRF}'_'${expo}'_genCorr.err

~/bolt \
--numThreads 20 \
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_scratch}'pheno_EC.txt \
--phenoCol '${expo}' \
--phenoCol '${CRF}' \
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
--statsFileBgenSnps '${expo_folder}'/'${CRF}'_'${expo}'_bgen.txt

' > ${expo_folder}/genCorr_${CRF}_${expo}.sh
sbatch ${expo_folder}/genCorr_${CRF}_${expo}.sh
echo Submitted genCorr for ${CRF} and ${expo}

done
done
done
