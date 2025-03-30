# Part of this project uses the summary stats from an external T2D GWAS
# However, that GWAS is in hg19, so we need to convert the rsIDs to hg38

dir_scratch='/n/groups/patel/nuno/PXS-pipeline/scratch/'
dir_T2DGWAS=${dir_scratch}'T2D/'
mkdir -p ${dir_T2DGWAS}
cd dir_T2DGWAS

# Installation of bigBedToBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bigBedToBed
chmod +x bigBedToBed
# downloading hg19 dbSNP data (contains rsID and hg19 coords) ~1.5Gb
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp153Common.bb
loc_snp='dbSnp153Common'
# use bigBedToBed to convert to bed file
./bigBedToBed ${loc_snp}.bb ${loc_snp}.bed
cut -f1-4 "${loc_snp}.bed" > "${loc_snp}_4col.bed"

# list of SNPs to extract
# for ease, I recommend using the ldsc-formatted sumstats
loc_LMM=${dir_scratch}'ldsc/behaviors/ldsc.BRS.sumstats.gz'
gunzip -c $loc_LMM | tail -n +2 | cut -f1 > rsIDs.txt
grep -Fwf rsIDs.txt ${loc_snp}_4col.bed > rsID_coordinates_hg19.bed # takes a while!



# T2D GWAS results have to be manually downloaded from the following site:
# https://www.diagram-consortium.org/downloads.html
# "T2DGGI GWAS all ancestry meta-analysis summary statistics"
# You must check the agreement box, download, and then specify zipped location here:
loc_T2DGWAS_zip='Suzuki.Nature2024.T2DGGI.MultiAncestry.sumstats'
unzip ${loc_T2DGWAS_zip}
loc_T2DGWAS='All_Metal_LDSC-CORR_Neff.v2.txt'

# the files that matter from this script are going to be:
# - the unzipped T2DGWAS file ('All_Metal_LDSC-CORR_Neff.v2.txt')
# - the resulting bed file: ('rsID_coordinates_hg19.bed')