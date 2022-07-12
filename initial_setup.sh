#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=8G
#SBATCH -o initial_setup.out
#SBATCH -e initial_setup.err

# Makes directory for each phenotype and 

loc_phenos="phenotypes.txt"
dir_out="/home/nur479/scratch3/PXS_pipeline/" #include last forward slash /
dir_script=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Creates ukkb_pheno and fields file that includes all phenotypes
#module load gcc/9.2.0 R/4.1.2
#R CMD BATCH initial_setup.R

lines=$(cat $loc_phenos)

#mkdir ${dir_out}results
for disease in $lines
do

cd ${dir_script}

subfolder=$(echo ${dir_out}${disease})
echo ${subfolder}
mkdir -p ${subfolder}
cp ${dir_script}/compute_PXS_LM_RF.R ${subfolder}/${disease}_compute_PXS_LM_RF.R

cd ${subfolder}

echo '#!/bin/sh
#SBATCH -c 1
#SBATCH -t 4-00:00
#SBATCH -p medium
#SBATCH --mem=32G
#SBATCH -o '${disease}'_compute_PXS_LM_RF.out
#SBATCH -e '${disease}'_compute_PXS_LM_RF.err

module load gcc/9.2.0 R/4.1.2
Rscript '${subfolder}'/'${disease}'_compute_PXS_LM_RF.R '${disease}'
' > ${subfolder}/${disease}_compute_PXS_LM_RF.sh
sbatch ${subfolder}/${disease}_compute_PXS_LM_RF.sh

done

