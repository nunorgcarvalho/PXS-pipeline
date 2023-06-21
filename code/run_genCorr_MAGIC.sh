dir_repo=${HOME}'/group_nuno/PXS-pipeline/'
dir_scratch=${dir_repo}'scratch/'
dir_MAGIC=${dir_scratch}'MAGIC/'
dir_script=${dir_repo}'code/'
dir_LD=${dir_scratch}LDscore/
dir_ldsc=${HOME}'/group_nuno/ldsc/'

## Downloading data, skip if already done
SKIP_DOWNLOADS='TRUE'

## Set which analyses to run
run_PXS="FALSE"
run_expos="TRUE"

if [ "$SKIP_DOWNLOADS" == FALSE ]; then
  cd ${dir_MAGIC}
  # downloads files
  wget https://magicinvestigators.org/downloads/MAGIC1000G_2hGlu_EUR.tsv.gz
  wget https://magicinvestigators.org/downloads/MAGIC1000G_FG_EUR.tsv.gz
  wget https://magicinvestigators.org/downloads/MAGIC1000G_FI_EUR.tsv.gz
  wget https://magicinvestigators.org/downloads/MAGIC1000G_HbA1c_EUR.tsv.gz
  # unzips them
  gunzip *.gz
  # renames them
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

if [ "$run_PXS" == TRUE ]; then # only runs if checked
  ## Converts PXS summary statistics to sumstats format
  loc_PXS_LMM=${dir_scratch}'T2D/LMM_T2D_bgen.txt'
  out_PXS="${loc_PXS_LMM%????}"
  sf_PXS=${out_PXS}.sumstats.gz
  
  if [ ! -f "$sf_PXS" ]; then # only runs if file does not already exist
    N=$(wc -l ${dir_scratch}'T2D/PXS_T2D.txt' | awk '{print $1}')
    echo ${loc_PXS_LMM}
    
    ${dir_ldsc}munge_sumstats.py \
    --sumstats ${loc_PXS_LMM} \
    --N ${N} \
    --snp SNP \
    --a1 ALLELE1 \
    --a2 ALLELE0 \
    --signed-sumstats BETA,0 \
    --p P_BOLT_LMM_INF \
    --out ${out_PXS}
    
    gunzip -c ${out_PXS}.sumstats.gz > ${out_PXS}.temp.tsv
    awk -F '\t' '{ print $1 "\t" $2 "\t" $3 }' ${out_PXS}.temp.tsv > ${out_PXS}.snplist
    rm ${out_PXS}.temp.tsv
  fi
fi

## Loops through the MAGIC phenotypes
MAGIC_sfs=$(ls -d ${dir_MAGIC}MAGIC_*.tsv)
sfs=($MAGIC_sfs)
length=${#sfs[@]}
for ((i=0; i<$length; i++)); do
  
  sf=${sfs[i]}
  out="${sf%????}"
  sf_munge=${out}.sumstats.gz
  echo 'Working on' ${sf_munge}
  
  if [ ! -f "$sf_munge" ]; then # only runs if file not already created
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
  fi
  
  if [ "$run_PXS" == TRUE ]; then # only runs if checked
    echo 'Running genetic correlation between' ${sf_munge} 'and PXS_T2D'
    # LD Score Regression
    ${dir_ldsc}ldsc.py \
    --rg ${sf_munge},${sf_PXS} \
    --ref-ld-chr ${dir_LD}LDscore. \
    --w-ld-chr ${dir_LD}LDscore. \
    --out ${out}_PXS_rg
  fi
  
  if [ "$run_expos" == TRUE ]; then # only runs if checked
    # loops through each exposure
    loc_exposures=$(echo ${dir_script}../input_data/exposures.txt)
    lines=$(cat $loc_exposures)
    echo 'All exposures:' ${lines[@]}
    for exposure in $lines; do
      echo 'Running loop for' ${exposure}
      dir_expo=${dir_scratch}exposures/${exposure}/
      loc_LMM_out=${dir_expo}${exposure}_BOLTLMM.out
      loc_sf_raw=${dir_expo}LMM_${exposure}_bgen.txt
      loc_out_sf="${loc_sf_raw%????}"
      loc_sf=${loc_out_sf}.sumstats.gz
      
      if [ ! -f "$loc_sf" ]; then
        # extracts sample size of GWAS
        line=$(grep "Number of individuals used in analysis: Nused =" $loc_LMM_out | head -n 1)
        N=$(echo $line | awk -F'= ' '{print $2}')
        # uses munge_stats.py from ldsc to properly format summary file for ldsc
        echo 'Formatting' ${exposure} 'summary file for ldsc'
        ${dir_ldsc}munge_sumstats.py \
        --sumstats ${loc_sf_raw} \
        --N ${N} \
        --snp SNP \
        --a1 ALLELE1 \
        --a2 ALLELE0 \
        --signed-sumstats BETA,0 \
        --p P_BOLT_LMM_INF \
        --out ${loc_out_sf}
      fi
      
      echo 'Running genetic correlation between' ${sf_munge} 'and' ${exposure}
      # LD Score Regression
      ${dir_ldsc}ldsc.py \
      --rg ${sf_munge},${loc_sf} \
      --ref-ld-chr ${dir_LD}LDscore. \
      --w-ld-chr ${dir_LD}LDscore. \
      --out ${out}_${exposure}_rg
    done
  fi
done