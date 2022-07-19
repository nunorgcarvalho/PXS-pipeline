# PXS Pipeline

## Order of scripts to run/edit:
1. `phenotypes.txt`: on each line ,write the name of the phenotypes (according to the coeffs and tally .csv) that an analysis should be performed on
2. `CRFs.txt`: table of clinical risk factors for each phenotype
3. `initial_setup.sh`: extracts necessary fields from phenotype file, and initiates scripts that compute the PXS for each phenotype
4. `run_BOLT_PXS.sh`: runs BOLT-LMM and BOLT-REML on the PXSs of the phenotypes
5. `run_BOLT_exposures.sh`: runs BOLT-LMM and BOLT-REML on the exposures of the PXS