#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-04:00
#SBATCH -p short
#SBATCH --mem=8G
#SBATCH -o geospatial_autocorrelations.out
#SBATCH -e geospatial_autocorrelations.err

# Makes list of exposures to analyze, appends PXSs to pheno file, and makes list of CRFs for each
module load gcc/9.2.0
module load R/4.1.2
module load gdal/3.1.4
module load geos/3.10.2
module load udunits/2.2.28

export UDUNITS2_INCLUDE=/n/app/udunits/2.2.28-gcc-9.2.0/include
export UDUNITS2_LIBS=/n/app/udunits/2.2.28-gcc-9.2.0/lib

R CMD BATCH geospatial_autocorrelations.R
