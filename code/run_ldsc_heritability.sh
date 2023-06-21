dir_repo=${HOME}'/group_nuno/PXS-pipeline/'
dir_scratch=${dir_repo}'scratch/'
dir_script=${dir_repo}'code/'
dir_LD=${dir_scratch}LDscore/
dir_ldsc=${HOME}'/group_nuno/ldsc/'

# loops through each exposure
loc_exposures=$(echo ${dir_script}../input_data/exposures.txt)
lines=$(cat $loc_exposures)
echo 'All exposures:' ${lines[@]}
for exposure in $lines
do

echo 'Running loop for' ${exposure}
dir_expo=${dir_scratch}exposures/${exposure}/
# extracts sample size of GWAS
loc_LMM_out=${dir_expo}${exposure}_BOLTLMM.out
line=$(grep "Number of individuals used in analysis: Nused =" $loc_LMM_out | head -n 1)
N=$(echo $line | awk -F'= ' '{print $2}')
# uses munge_stats.py from ldsc to properly format summary file for ldsc
loc_sf_raw=${dir_expo}LMM_${exposure}_bgen.txt
loc_out_sf="${loc_sf_raw%????}"
loc_sf=${loc_out_sf}.sumstats.gz

if [ ! -f "$loc_sf" ]; then
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
fi

loc_out=${dir_expo}ldsc_h2_${exposure}
# LD Score Regression. Quick (< 1 minute)
echo 'Running ldsc h2 estimation'
${dir_ldsc}ldsc.py \
--h2 ${loc_sf} \
--ref-ld-chr ${dir_LD}LDscore. \
--w-ld-chr ${dir_LD}LDscore. \
--out ${loc_out}
done