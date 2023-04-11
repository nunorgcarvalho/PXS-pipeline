dir_MAGIC='~/scratch3/key_data/MAGIC'
dir_script='~/jobs/PXS_pipeline/code'
dir_scratch='/home/nur479/group_nuno/PXS_pipeline/scratch/'

## Downloading data, skip if already done
SKIP_DOWNLOADS='true'
if [ "$SKIP_DOWNLOADS" == false ]; then
cd ${dir_MAGIC}

wget https://magicinvestigators.org/downloads/PROIadjBMI_META_GWAS_CAT_20221113.tsv.gz
wget https://magicinvestigators.org/downloads/PROInoBMI_META_GWAS_CAT_20221113.tsv.gz
wget https://magicinvestigators.org/downloads/MAGIC1000G_2hGlu_EUR.tsv.gz
wget https://magicinvestigators.org/downloads/MAGIC1000G_FG_EUR.tsv.gz
wget https://magicinvestigators.org/downloads/MAGIC1000G_FI_EUR.tsv.gz
wget https://magicinvestigators.org/downloads/MAGIC1000G_HbA1c_EUR.tsv.gz

gunzip *.gz

cd ${dir_script}
fi

loc_PXS_LMM=${}