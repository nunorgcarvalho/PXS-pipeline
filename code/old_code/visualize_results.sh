#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-04:00
#SBATCH -p short
#SBATCH --mem=12G
#SBATCH -o visualize_results.out
#SBATCH -e visualize_results.err

# Run this after running:
# initial_setup.sh
# run_BOLT_PXS.sh
# run_BOLT_exposures.sh
# run_BOLT_genCorr_CRFs.sh
# aggregate_results.sh

dir_script="/home/nur479/jobs/PXS_pipeline/code/"
cd ${dir_script}

module load gcc/9.2.0 R/4.1.2
R CMD BATCH visualize_results.R