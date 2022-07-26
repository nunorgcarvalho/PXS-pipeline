#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-04:00
#SBATCH -p short
#SBATCH --mem=12G
#SBATCH -o aggregate_results.out
#SBATCH -e aggregate_results.err

dir_script="/home/nur479/jobs/PXS_pipeline/"
cd ${dir_script}

module load gcc/9.2.0 R/4.1.2
R CMD BATCH aggregate_results.R
