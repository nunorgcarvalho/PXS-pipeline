dir_repo=${HOME}'/group_nuno/PXS-pipeline/'
dir_scratch=${dir_repo}'scratch/'
dir_MAGIC=${dir_scratch}'MAGIC/'
dir_script=${dir_repo}'code/'
dir_LD=${dir_scratch}LDscore/
dir_ldsc=${HOME}'/group_nuno/ldsc/'

## Downloading data, skip if already done
SKIP_DOWNLOADS='true'
if [ "$SKIP_DOWNLOADS" == false ]; then
cd ${dir_MAGIC}
# downloads files
#wget https://magicinvestigators.org/downloads/PROIadjBMI_META_GWAS_CAT_20221113.tsv.gz
#wget https://magicinvestigators.org/downloads/PROInoBMI_META_GWAS_CAT_20221113.tsv.gz
wget https://magicinvestigators.org/downloads/MAGIC1000G_2hGlu_EUR.tsv.gz
wget https://magicinvestigators.org/downloads/MAGIC1000G_FG_EUR.tsv.gz
wget https://magicinvestigators.org/downloads/MAGIC1000G_FI_EUR.tsv.gz
wget https://magicinvestigators.org/downloads/MAGIC1000G_HbA1c_EUR.tsv.gz
# unzips them
gunzip *.gz
# renames them
#mv PROIadjBMI_META_GWAS_CAT_20221113.tsv MAGIC_FPIadj.tsv
#mv PROInoBMI_META_GWAS_CAT_20221113.tsv MAGIC_FPInoadj.tsv
mv MAGIC1000G_2hGlu_EUR.tsv MAGIC_2hGlu.tsv
mv MAGIC1000G_FG_EUR.tsv MAGIC_FG.tsv
mv MAGIC1000G_FI_EUR.tsv MAGIC_FI.tsv
mv MAGIC1000G_HbA1c_EUR.tsv MAGIC_HbA1c.tsv

# Downloads LD scores for Europeans
cd ${dir_scratch}
wget https://zenodo.org/record/7796478/files/1000G_Phase3_ldscores.tgz
tar -zxvf 1000G_Phase3_ldscores.tgz
rm 1000G_Phase3_ldscores.tgz

cd ${dir_script}
fi

## Converts PXS summary statistics to sumstats format
loc_PXS_LMM=${dir_scratch}'T2D/LMM_T2D_bgen.txt'
N=$(wc -l ${dir_scratch}'T2D/PXS_T2D.txt' | awk '{print $1}')
#N=282118
echo ${loc_PXS_LMM}
out_PXS="${loc_PXS_LMM%????}"

${dir_ldsc}munge_sumstats.py \
--sumstats ${loc_PXS_LMM} \
--N ${N} \
--snp SNP \
--a1 ALLELE1 \
--a2 ALLELE0 \
--signed-sumstats BETA,0 \
--p P_BOLT_LMM_INF \
--out ${out_PXS}

sf_PXS=${out_PXS}.sumstats
gunzip -c ${out_PXS}.sumstats.gz > ${out_PXS}.temp.tsv
awk -F '\t' '{ print $1 "\t" $2 "\t" $3 }' ${out_PXS}.temp.tsv > ${out_PXS}.snplist
rm ${out_PXS}.temp.tsv
#gunzip ${sf_PXS}.gz


## Loops through the MAGIC phenotypes
MAGIC_sfs=$(ls -d ${dir_MAGIC}MAGIC_*.tsv)
sfs=($MAGIC_sfs)
length=${#sfs[@]}
for ((i=0; i<$length; i++)); do

sf=${sfs[i]}
out="${sf%????}"

echo $sf
# Converts MAGIC summary statistics to sumstats format
${dir_ldsc}munge_sumstats.py \
--sumstats $sf \
--N-col sample_size \
--snp variant \
--a1 effect_allele \
--a2 other_allele \
--merge-alleles ${out_PXS}.snplist \
--chunksize 500000 \
--signed-sumstats beta,0 \
--out ${out}

sf_munge=${out}.sumstats
# LD Score Regression
${dir_ldsc}ldsc.py \
--rg ${out}.sumstats.gz,${out_PXS}.sumstats.gz \
--ref-ld-chr ${dir_LD}LDscore. \
--w-ld-chr ${dir_LD}LDscore. \
--out ${out}_PXS_rg


done
