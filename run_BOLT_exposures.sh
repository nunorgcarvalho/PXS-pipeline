# sets directories and paths
loc_exposures="exposures.txt" # will be created by make_exposures_list.R by this script
dir_out="/home/nur479/scratch3/PXS_pipeline/" #include last forward slash /
dir_script="/home/nur479/jobs/PXS_pipeline"

# Makes list of exposures to analyze
# Slow to run from non-interactive node, but it is doable (about 2 minutes)
module load gcc/9.2.0 R/4.1.2
R CMD BATCH make_exposures_list.R
echo 'Made list of exposures'

dir_exposures=$(echo ${dir_out}exposures)
mkdir -p ${dir_exposures}

lines=$(cat $loc_exposures)
for exposure in $lines
do

subfolder=$(echo ${dir_exposures}/${exposure})
mkdir -p ${subfolder}
cd ${subfolder}

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
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_out}'bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_out}'pheno_EC.txt \
--phenoCol '${exposure}' \
--covarFile '${dir_out}'pheno_EC.txt \
--covarCol f.31.0.0 \
--covarCol f.54.0.0 \
--qCovarCol f.34.0.0 \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--lmm \
--verboseStats \
--statsFile '${subfolder}'/LMM_'${exposure}'.txt \
--bgenFile /n/no_backup2/patel/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${subfolder}'/LMM_'${exposure}'_bgen.txt

' > ${subfolder}/${exposure}_BOLTLMM.sh
sbatch ${subfolder}/${exposure}_BOLTLMM.sh
echo 'Submitted BOLT-LMM for '${exposure}

#########################################
## Creates BOLT-REML script and submits##
#########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o '${exposure}'_BOLTREML.out
#SBATCH -e '${exposure}'_BOLTREML.err

~/bolt \
--numThreads 20 \
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_out}'bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_out}'pheno_EC.txt \
--phenoCol '${exposure}' \
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
--statsFileBgenSnps '${subfolder}'/LMM_'${exposure}'_bgen.txt

' > ${subfolder}/${exposure}_BOLTREML.sh
sbatch ${subfolder}/${exposure}_BOLTREML.sh
echo 'Submitted BOLT-REML for '${exposure}

done