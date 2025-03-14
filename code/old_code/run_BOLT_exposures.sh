# sets directories and paths (include last forward slash / )
dir_script="/home/nur479/group_nuno/PXS-pipeline/code/"
dir_scratch="/home/nur479/group_nuno/PXS-pipeline/scratch/"

dir_exposures=$(echo ${dir_scratch}exposures)
mkdir -p ${dir_exposures}

loc_exposures=$(echo ${dir_script}../input_data/exposures.txt)
lines=$(cat $loc_exposures)
for exposure in $lines
do

subfolder=$(echo ${dir_exposures}/${exposure})
mkdir -p ${subfolder}
cd ${subfolder}

# skips if already calculated
out_file=${subfolder}/${exposure}_BOLTLMM.out
if [[ -f "$out_file" ]]; then
    continue
fi

########################################
## Creates BOLT-LMM script and submits##
########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-00:00
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o '${exposure}'_BOLTLMM.out
#SBATCH -e '${exposure}'_BOLTLMM.err

~/bolt \
--numThreads 20 \
--bed /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/no_backup2/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_scratch}'pheno_EC.txt \
--phenoCol '${exposure}' \
--covarFile '${dir_scratch}'pheno_EC.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--lmm \
--verboseStats \
--statsFile '${subfolder}'/LMM_'${exposure}'.txt \
--bgenFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/uk_biobank/ukb_genetics/22881/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/LMM_'${exposure}'_bgen.txt

' > ${subfolder}/${exposure}_BOLTLMM.sh
sbatch ${subfolder}/${exposure}_BOLTLMM.sh
echo 'Submitted BOLT-LMM for '${exposure}

done
