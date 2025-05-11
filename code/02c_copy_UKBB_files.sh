# copies UKBB genotype files onto my scratch so I can do analyses on it
dir_source='/n/no_backup2/patel/uk_biobank/main_data_9512/'
dir_target='/n/scratch/users/n/nur479/UKBB_geno_files/'
file_prefix='UKBB_geno_chr'


for chr in {1..22}; do
    echo "Copying chromosome ${chr} bed/bim/fam"
    cp -v ${dir_source}ukb_cal_chr${chr}_v2.bed ${dir_target}${file_prefix}${chr}.bed
    cp ${dir_source}ukb_snp_chr${chr}_v2.bim ${dir_target}${file_prefix}${chr}.bim
    cp ${dir_source}ukb2288_cal_chr${chr}_v2_s488366.fam ${dir_target}${file_prefix}${chr}.fam
done
