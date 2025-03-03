dir_repo=${HOME}'/group_nuno/PXS-pipeline/'
dir_scratch=${dir_repo}'scratch/'
dir_script=${dir_repo}'code/'
dir_LD=${dir_scratch}LDscore/
dir_ldsc=${HOME}'/group_nuno/ldsc/'
dir_T2D=${dir_scratch}'T2D/'
# set which software to use for gencorr
run_ldsc="FALSE"
run_BOLTREML="TRUE"

# LDsc gencorr
if [ "$run_ldsc" == "TRUE" ]; then
# corrects formatting for all the two T2Ds (PXS_T2D already done in run_genCorr_MAGIC.sh)
for pheno in T2D_onset T2D_all
do
echo 'Running loop for' ${pheno}

loc_LMM_out=${dir_T2D}${pheno}_BOLTLMM.out
line=$(grep "Number of individuals used in analysis: Nused =" $loc_LMM_out | head -n 1)
N=$(echo $line | awk -F'= ' '{print $2}')
# uses munge_stats.py from ldsc to properly format summary file for ldsc
loc_sf_raw=${dir_T2D}LMM_${pheno}_bgen.txt
loc_out_sf="${loc_sf_raw%????}"
# takes about 3 minutes
echo 'Formatting summary file for ldsc'
${dir_ldsc}munge_sumstats.py \
--sumstats ${loc_sf_raw} \
--N ${N} \
--snp SNP \
--a1 ALLELE1 \
--a2 ALLELE0 \
--signed-sumstats BETA,0 \
--p P_BOLT_LMM_INF \
--out ${loc_out_sf}

done

phenos1=(PXS_T2D PXS_T2D T2D_onset)
phenos2=(T2D_onset T2D_all T2D_all)
# computes gencorr for 3 pairs of PXS_T2D, T2D_onset, T2D_all
for i in 0 1 2
do
# gets names and paths to formatted summary files
pheno1=${phenos1[i]}
pheno2=${phenos2[i]}

pheno1_sf=${dir_T2D}LMM_${pheno1}_bgen.sumstats.gz
pheno2_sf=${dir_T2D}LMM_${pheno2}_bgen.sumstats.gz

# runs ldsc rg estimate
echo 'Running gencorr for' ${pheno1} 'x' ${pheno2}
${dir_ldsc}ldsc.py \
--rg ${pheno1_sf},${pheno2_sf} \
--ref-ld-chr ${dir_LD}LDscore. \
--w-ld-chr ${dir_LD}LDscore. \
--out ${dir_T2D}LDsc_gencorr_${pheno1}_${pheno2}

done

fi

# BOLT-REML gencorr
if [ "$run_BOLTREML" == "TRUE" ]; then

cd ${dir_T2D}

# genCorr(PXS_T2D, T2D_onset)

#########################################
## Creates BOLT-REML script and submits##
#########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 2-11:59
#SBATCH -p medium
#SBATCH --mem=100G
#SBATCH -o T2D_PXS_onset_genCorr.out
#SBATCH -e T2D_PXS_onset_genCorr.err

~/bolt \
--numThreads 20 \
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_scratch}'pheno_EC.txt \
--phenoCol PXS_T2D \
--phenoCol T2D_onset \
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
--statsFileBgenSnps '${dir_T2D}'_PXS_onset_bgen.txt

' > ${dir_T2D}/genCorr_T2D_PXS_onset.sh
#sbatch ${dir_T2D}/genCorr_T2D_PXS_onset.sh
echo Submitted genCorr for PXS_T2D and T2D_onset


# genCorr(PXS_T2D, PXS_T2D_BMIadj)

#########################################
## Creates BOLT-REML script and submits##
#########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-11:59
#SBATCH -p medium
#SBATCH --mem=100G
#SBATCH -o T2D_PXS_BMIadj2_genCorr.out
#SBATCH -e T2D_PXS_BMIadj2_genCorr.err

~/bolt \
--numThreads 20 \
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_script}'../input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_scratch}'pheno_EC.txt \
--phenoCol PXS_T2D \
--phenoCol PXS_T2D_BMIadj \
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
--statsFileBgenSnps '${dir_T2D}'_PXS_BMIadj2_bgen.txt

' > ${dir_T2D}/genCorr_T2D_PXS_BMIadj2.sh
sbatch ${dir_T2D}/genCorr_T2D_PXS_BMIadj2.sh
echo Submitted genCorr for PXS_T2D and PXS_T2D_BMIadj2

fi
