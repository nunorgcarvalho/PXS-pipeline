# PXS Pipeline

## Input files to edit:
- `phenotypes.txt`: on each line, write the name of the phenotypes (according to the coeffs and tally .csv) that an analysis should be performed on
- `CRFs_table.txt`: table of clinical risk factors for each phenotype

## Order of scripts to run:
1. `initial_setup.sh`: extracts necessary fields from phenotype file, and initiates scripts that compute the PXS for each phenotype
2. `run_BOLT_PXS.sh`: runs BOLT-LMM and BOLT-REML on the PXSs of the phenotypes
3. `run_BOLT_exposures.sh`: runs BOLT-LMM and BOLT-REML on the exposures of the PXS
4. `run_genCorr_CRFs.sh`: run BOLT-REML (genetic correlation) for each PXS and its disease's clinical risk factors

## Visualizing/Summarizing data:
- `aggregate_results.R`: extracts important information from BOLT-REML results for XWAS diseases and exposures
- `geographical_clustering.R`: calculates ANOVA and saves plots for PXS and exposures per assessment center