#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=16G
#SBATCH -o initial_setup.out
#SBATCH -e initial_setup.err


# sets directories and paths (include last forward slash / )
dir_script="/home/nur479/group_nuno/PXS-pipeline/code/"
dir_scratch="/home/nur479/group_nuno/PXS-pipeline/scratch/"

cd ${dir_script}

# Creates ukkb_pheno and fields file that includes all phenotypes
module load gcc/9.2.0 R/4.1.2
R CMD BATCH generate_PXS_coeffs.R

loc_phenolist=$(echo ${dir_script}../input_data/phenotypes.txt)
phenos=$(cat $loc_phenolist)

for disease in $phenos
do

cd ${dir_script}

subfolder=$(echo ${dir_scratch}${disease})
echo ${subfolder}
mkdir -p ${subfolder}
cp ${dir_script}compute_PXS_LM.R ${subfolder}/${disease}_compute_PXS_LM.R

#cd ${subfolder}

echo '#!/bin/sh
#SBATCH -c 1
#SBATCH -t 0-02:00
#SBATCH -p short
#SBATCH --mem=32G
#SBATCH -o '${disease}'_compute_PXS_LM.out
#SBATCH -e '${disease}'_compute_PXS_LM.err

module load gcc/9.2.0 R/4.1.2
Rscript '${subfolder}'/'${disease}'_compute_PXS_LM.R '${disease}'
' > ${subfolder}/${disease}_compute_PXS_LM.sh
sbatch ${subfolder}/${disease}_compute_PXS_LM.sh

done

